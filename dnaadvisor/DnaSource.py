from copy import copy
from .optimization import (optimize_cuts_with_graph_twostep,
                           NoSolutionFoundError)
from dnachisel import Constraint, DnaCanvas
from DnaOrderingPlan import DnaQuote, DnaOrderingPlan

class DnaSource:

    def get_quote(self, sequence, max_delay=None, max_price=None,
                  with_ordering_plan=False, time_resolution=1.0):

        if self.memoize:
            args = (sequence, max_delay, max_price, with_ordering_plan)
            quote = self.memoize_dict.get(args, None)
            if quote is not None:
                return quote

        if not self.verify_constraints(sequence):
            quote = DnaQuote(self, sequence, accepted=False,
                             message="Constraints do not pass")

        elif max_price is not None:
            quote = self.get_best_delay_under_price_limit(
                sequence,
                max_price=max_price,
                time_resolution=time_resolution,
                with_ordering_plan=with_ordering_plan,

            )
        else:
            quote = self.get_best_price(
                sequence,
                max_delay=max_delay,
                with_ordering_plan=with_ordering_plan
            )

        if self.memoize:
            self.memoize_dict[args] = quote

        return quote

    def __repr__(self):
        return self.name

    def get_best_delay_under_price_limit(self, sequence, max_price,
                                         time_resolution,
                                         with_ordering_plan=None):
        """dichotomize, kind of"""
        def f(max_delay):
            return self.get_quote(sequence, max_delay=max_delay,
                                  with_ordering_plan=with_ordering_plan)
        quote = f(None)
        if (not quote.accepted) or (quote.price > max_price):
            return DnaQuote(
                self, sequence, accepted=False,
                message="Price over limit even without time limit."
            )
        step_size = quote.delay / 2.0
        delay = quote.delay - step_size
        best_quote = quote
        while step_size > time_resolution:
            quote = f(delay)
            if quote.accepted and (quote.price <= max_price):
                best_quote = quote
                delay -= step_size
            else:
                delay += step_size
            step_size /= 2.0

        return best_quote


    def verify_constraints(self, sequence):
        constraints = self.sequence_constraints
        dnachisel_constraints = [
            constraint for constraint in constraints
            if isinstance(constraint, Constraint)
        ]

        if dnachisel_constraints != []:
            canvas = DnaCanvas(sequence, dnachisel_constraints)
            constraints = [
                constraint for constraint in constraints
                if not isinstance(constraint, Constraint)
            ] + [lambda seq: canvas.all_constraints_pass()]

        return all(constraint(sequence) for constraint in constraints)




class DnaSourcesComparator(DnaSource):

    def __init__(self, dna_sources, memoize=False):
        self.dna_sources = dna_sources
        self.sequence_constraints = ()
        self.memoize = memoize
        self.memoize_dict = {}

    def get_best_price(self, sequence, max_delay=None,
                       with_ordering_plan=False):
        quotes = [
            source.get_quote(sequence, max_delay=max_delay,
                             with_ordering_plan=with_ordering_plan)
            for source in self.dna_sources
        ]
        accepted_quotes = [quote for quote in quotes if quote.accepted]
        if accepted_quotes == []:
            return DnaQuote(self, sequence, accepted=False,
                            message="Found no source accepting the sequence.")
        else:
            return min(accepted_quotes, key=lambda quote: quote.price)


class DnaAssemblyStation(DnaSource):

    def __init__(self, name, dna_assembly_method, dna_source, memoize=False,
                 **solve_kwargs):
        self.name = name
        self.dna_assembly_method = dna_assembly_method
        self.dna_source = dna_source
        self.solve_kwargs = solve_kwargs
        self.extra_delay = dna_assembly_method.duration
        self.extra_price = dna_assembly_method.cost
        self.sequence_constraints = dna_assembly_method.sequence_constraints
        self.memoize = memoize
        self.memoize_dict = {}

    def get_quote_for_sequence_segment(self, sequence, segment, max_delay=None,
                                       **kwargs):
        fragment_to_order = self.dna_assembly_method.compute_fragment_sequence(
            sequence, segment, **kwargs
        )
        return self.dna_source.get_quote(fragment_to_order,
                                         max_delay=max_delay)

    def get_ordering_plan_from_cuts(self, sequence, cuts, max_delay=None):
        cuts = sorted(cuts)
        return DnaOrderingPlan({
            segment: self.get_quote_for_sequence_segment(sequence, segment,
                                                         max_delay=max_delay)
            for segment in zip(cuts, cuts[1:])
        })

    def get_ordering_plan_for_sequence(self, sequence, max_delay=None,
                                       return_graph=False,
                                       nucleotide_resolution=None,
                                       refine_resolution=None,
                                       progress_bars=True):
        def segment_score(segment):
            quote = self.get_quote_for_sequence_segment(
                sequence, segment, max_delay=max_delay
            )
            if not quote.accepted:
                return -1
            else:
                return quote.price
        assembly = self.dna_assembly_method
        if nucleotide_resolution is None:
            nucleotide_resolution = \
                self.solve_kwargs.get("nucleotide_resolution", None)
        if refine_resolution is None:
            refine_resolution = \
                self.solve_kwargs.get("refine_resolution", None)
        graph, best_cuts = optimize_cuts_with_graph_twostep(
            sequence_length=len(sequence),
            segment_score_function=segment_score,
            location_filters=assembly.location_filters,
            segment_filters=assembly.segment_filters,
            min_segment_length=assembly.min_segment_length,
            max_segment_length=assembly.max_segment_length,
            forced_cuts=assembly.force_cuts(sequence),
            initial_resolution=nucleotide_resolution,
            refine_resolution=refine_resolution,
            progress_bars=progress_bars
        )
        ordering_plan = self.get_ordering_plan_from_cuts(sequence, best_cuts,
                                                         max_delay=max_delay)
        if return_graph:
            return graph, ordering_plan
        else:
            return ordering_plan

    def get_best_price(self, sequence, max_delay=None,
                       with_ordering_plan=False):

        if (max_delay is not None):
            max_delay = max_delay - self.extra_delay
            if max_delay < 0:
                return DnaQuote(self, sequence, accepted=False,
                                message="Delay too short")
        try:
            ordering_plan = self.get_ordering_plan_for_sequence(
                sequence, max_delay=max_delay, **self.solve_kwargs
            )
        except NoSolutionFoundError:
            return DnaQuote(self, sequence, accepted=False,
                            message="No solution found !")
        total_price = ordering_plan.total_price() + self.extra_price
        delay = ordering_plan.overall_delay()
        total_delay = None if delay is None else delay + self.extra_delay
        if not with_ordering_plan:
            ordering_plan = None
        return DnaQuote(self, sequence, price=total_price, delay=total_delay,
                        ordering_plan=ordering_plan)


class ExternalDnaOffer(DnaSource):

    def __init__(self, name, sequence_constraints, price_function, delay=None,
                 memoize=False):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.price_function = price_function
        self.delay_function = delay if callable(delay) else (lambda *a: delay)
        self.memoize = memoize
        self.memoize_dict = {}

    def __repr__(self):
        return self.name

    def get_best_price(self, sequence, max_delay=None,
                       with_ordering_plan=False):

        delay = self.delay_function(sequence)
        price = self.price_function(sequence)
        accepted = (max_delay is None) or (delay <= max_delay)
        return DnaQuote(self, sequence, accepted=accepted, delay=delay,
                        price=price)

    def get_best_delay_under_price_limit(self, sequence, max_price):

        delay = self.delay_function(sequence)
        price = self.price_function(sequence)
        return DnaQuote(self, sequence, price=price, delay=delay,
                        accepted=price < max_price)

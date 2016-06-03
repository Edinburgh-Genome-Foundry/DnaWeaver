from copy import copy
from .optimization import (optimize_cuts_with_graph_twostep,
                           NoSolutionFoundError)
from dnachisel import Constraint, DnaCanvas
from DnaOrderingPlan import DnaQuote, DnaOrderingPlan


class DnaSource:

    def get_quote(self, sequence, max_lead_time=None, max_price=None,
                  with_ordering_plan=False, time_resolution=1.0):

        if self.memoize:
            args = (sequence, max_lead_time, max_price, with_ordering_plan)
            quote = self.memoize_dict.get(args, None)
            if quote is not None:
                return quote

        if not self.verify_constraints(sequence):
            quote = DnaQuote(self, sequence, accepted=False,
                             message="Constraints do not pass")

        elif max_price is not None:
            quote = self.get_best_lead_time_under_price_limit(
                sequence,
                max_price=max_price,
                time_resolution=time_resolution,
                with_ordering_plan=with_ordering_plan,

            )
        else:
            quote = self.get_best_price(
                sequence,
                max_lead_time=max_lead_time,
                with_ordering_plan=with_ordering_plan
            )

        if self.memoize:
            self.memoize_dict[args] = quote

        return quote

    def __repr__(self):
        return self.name

    def get_best_lead_time_under_price_limit(self, sequence, max_price,
                                             time_resolution,
                                             with_ordering_plan=None):
        """dichotomize, kind of"""
        def f(max_lead_time):
            return self.get_quote(sequence, max_lead_time=max_lead_time,
                                  with_ordering_plan=with_ordering_plan)
        quote = f(None)
        if (not quote.accepted) or (quote.price > max_price):
            return DnaQuote(
                self, sequence, accepted=False,
                message="Price over limit even without time limit."
            )
        step_size = quote.lead_time / 2.0
        lead_time = quote.lead_time - step_size
        best_quote = quote
        while step_size > time_resolution:
            quote = f(lead_time)
            if quote.accepted and (quote.price <= max_price):
                best_quote = quote
                lead_time -= step_size
            else:
                lead_time += step_size
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

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):
        quotes = [
            source.get_quote(sequence, max_lead_time=max_lead_time,
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

    def __init__(self, name, assembly_method, dna_source, memoize=False,
                 **solve_kwargs):
        self.name = name
        self.assembly_method = assembly_method
        self.dna_source = dna_source
        self.solve_kwargs = solve_kwargs
        self.extra_time = assembly_method.duration
        self.extra_price = assembly_method.cost
        self.sequence_constraints = assembly_method.sequence_constraints
        self.memoize = memoize
        self.memoize_dict = {}

    def get_quote_for_sequence_segment(self, sequence, segment,
                                       max_lead_time=None,
                                       **kwargs):
        fragment_to_order = self.assembly_method.compute_fragment_sequence(
            sequence, segment, **kwargs
        )
        return self.dna_source.get_quote(fragment_to_order,
                                         max_lead_time=max_lead_time)

    def get_ordering_plan_from_cuts(self, sequence, cuts, max_lead_time=None):
        cuts = sorted(cuts)
        return DnaOrderingPlan({
            segment: self.get_quote_for_sequence_segment(
                sequence, segment,
                max_lead_time=max_lead_time
            )
            for segment in zip(cuts, cuts[1:])
        })

    def get_ordering_plan_for_sequence(self, sequence, max_lead_time=None,
                                       return_graph=False,
                                       nucleotide_resolution=None,
                                       refine_resolution=None,
                                       a_star_factor=0,
                                       progress_bars=True):
        def segment_score(segment):
            quote = self.get_quote_for_sequence_segment(
                sequence, segment, max_lead_time=max_lead_time
            )
            if not quote.accepted:
                return -1
            else:
                return quote.price
        assembly = self.assembly_method
        if nucleotide_resolution is None:
            nucleotide_resolution = \
                self.solve_kwargs.get("nucleotide_resolution", 1)
        if refine_resolution is None:
            refine_resolution = \
                self.solve_kwargs.get("refine_resolution", 1)

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
            progress_bars=progress_bars,
            a_star_factor=a_star_factor
        )
        ordering_plan = self.get_ordering_plan_from_cuts(
            sequence, best_cuts,
            max_lead_time=max_lead_time
        )
        if return_graph:
            return graph, ordering_plan
        else:
            return ordering_plan

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):

        if (max_lead_time is not None):
            max_lead_time = max_lead_time - self.extra_time
            if max_lead_time < 0:
                return DnaQuote(self, sequence, accepted=False,
                                message="Time limit too short")
        try:
            ordering_plan = self.get_ordering_plan_for_sequence(
                sequence, max_lead_time=max_lead_time, **self.solve_kwargs
            )
        except NoSolutionFoundError:
            return DnaQuote(self, sequence, accepted=False,
                            message="No solution found !")
        total_price = ordering_plan.total_price() + self.extra_price
        lead_time = ordering_plan.overall_lead_time()
        total_duration = (None if lead_time is None else
                          lead_time + self.extra_time)
        if not with_ordering_plan:
            ordering_plan = None
        return DnaQuote(self, sequence, price=total_price,
                        lead_time=total_duration,
                        ordering_plan=ordering_plan)


class ExternalDnaOffer(DnaSource):

    def __init__(self, name, sequence_constraints, price_function,
                 lead_time=None, memoize=False):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.price_function = price_function
        self.lead_time = (lead_time if callable(lead_time)
                          else (lambda *a: lead_time))
        self.memoize = memoize
        self.memoize_dict = {}

    def __repr__(self):
        return self.name

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):

        lead_time = self.lead_time(sequence)
        price = self.price_function(sequence)
        accepted = (max_lead_time is None) or (lead_time <= max_lead_time)
        return DnaQuote(self, sequence, accepted=accepted,
                        lead_time=lead_time, price=price)

    def get_best_lead_time_under_price_limit(self, sequence, max_price):

        lead_time = self.lead_time(sequence)
        price = self.price_function(sequence)
        return DnaQuote(self, sequence, price=price, lead_time=lead_time,
                        accepted=price < max_price)


class PcrOutStation(DnaSource):

    def __init__(self, name, primers_dna_source, blast_database=None,
                 sequences=None, pcr_homology_length=25,
                 max_overhang_length=40, extra_cost=0, extra_time=0,
                 max_amplicon_length=None, blast_word_size=50, memoize=False,
                 sequence_constraints=()):
        self.name = name
        self.blast_database = blast_database
        self.primers_dna_source = primers_dna_source
        self.pcr_homology_length = pcr_homology_length
        self.max_overhang_length = max_overhang_length
        self.extra_time = extra_time
        self.extra_cost = extra_cost
        self.max_amplicon_length = max_amplicon_length
        self.blast_word_size = blast_word_size
        if max_amplicon_length is not None:
            sequence_constraints = ([lambda seq: len(seq) <
                                     max_amplicon_length] + list(constraints))
        self.sequence_constraints = sequence_constraints
        self.memoize = memoize
        self.memoize_dict = {}
        self.sequences = sequences

    def get_hits(self, sequence):
        if self.sequences is not None:
            result = []
            for dna_name, seq in self.sequences:
                index = seq.find(sequence)
                if index != -1:
                    result.append((dna_name, (index, index+len(sequence))))
            return result
        else:
            record = blast_sequence(sequence, self.blast_database,
                                    perc_identity=100,
                                    word_size=self.blast_word_size)
            return [
                (hit.sbjct, (hit.query_start, hit.query_end))
                for alignment in record.alignments
                for hit in alignment.hsps
            ]


    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):

        for subject, (hit_start, hit_end) in self.get_hits(sequence):

            largest_overhang = max(hit_start, len(sequence) - hit_end)

            if largest_overhang <= self.max_overhang_length:
                primer_l_end = hit_start + self.pcr_homology_length
                primer_left = sequence[:primer_l_end]
                primer_r_end = hit_end - self.pcr_homology_length
                primer_right = reverse_complement(sequence[primer_r_end:])

                primer_max_lead_time = (None if max_lead_time is None else
                                     max_lead_time - self.extra_duration)
                quotes = [
                    self.primers_dna_source.get_quote(
                        primer, max_lead_time=primer_max_lead_time
                    )
                    for primer in [primer_left, primer_right]
                ]
                if not all(quote.accepted for quote in quotes):
                    continue  # primers inorderable

                if max_lead_time is not None:
                    overall_lead_time = (max(quote.lead_time
                                             for quote in quotes) +
                                         self.extra_time)
                else:
                    overall_lead_time = None
                total_price = (sum(quote.price for quote in quotes) +
                               self.extra_cost)

                if with_ordering_plan:
                    ordering_plan = DnaOrderingPlan({
                        (0, primer_left_end): quotes[0],
                        (primer_right_end, len(sequence)): quotes[0]
                    })
                else:
                    ordering_plan = None
                return DnaQuote(self, sequence, accepted=True,
                                lead_time=overall_lead_time,
                                price=total_price,
                                ordering_plan=ordering_plan)

        return DnaQuote(self, sequence, accepted=False,
                        message="No valid match found")

    def pre_blast(self, sequence):
        self.sequences = None
        self.sequences = [
            (subject, sequence[start:end])
            for subject, (start, end) in self.get_hits(sequence)
        ]
        print len(self.sequences)

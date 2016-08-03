from copy import copy
from .optimization import (optimize_cuts_with_graph_twostep,
                           NoSolutionFoundError)

from biotools import blast_sequence, largest_common_substring, reverse_complement
from DnaOrderingPlan import DnaQuote, DnaOrderingPlan
import numpy as np
import networkx as nx
from collections import defaultdict

# Attempt import of optional module DNA Chisel
try:
    from dnachisel import Constraint, DnaCanvas
    DNACHISEL_AVAILABLE = True
except:
    DNACHISEL_AVAILABLE = False


class DnaSource:
    """Base class for all DnaSources, which are the elements of the supply
    networks used to define assembly problems in DnaWeaver."""

    def get_quote(self, sequence, max_lead_time=None, max_price=None,
                  with_ordering_plan=False, time_resolution=1.0):
        """Return a DnaQuote with price, lead time, etc. for a given sequence.

        Parameters
        ----------


        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times
          Note that either this parameter or `max_lead_time` must be None

        with_ordering_plan
          If True, the ordering plan is added to the quote

        time_resolution
          Time resolution for the bisecting search if `max_price` is not None.

        Returns
        -------

        A DnaQuote object.

        """

        if max_lead_time is None:
            max_lead_time = np.inf

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
        """Return the quote with fastest lead time under the budget constraint

        Parameters
        ----------



        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times
          Note that either this parameter or `max_lead_time` must be None

        with_ordering_plan
          If True, the ordering plan is added to the quote

        time_resolution
          Time resolution for the bisecting search if `max_price` is not None.

        """
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
        """Return True iff `sequence` passes all `self.sequence_constraints`

        Will automatically process DNA-Chisel constraints that would be in
        `self.sequence_constraints`

        """
        constraints = self.sequence_constraints
        dnachisel_constraints = [
            constraint for constraint in constraints
            if isinstance(constraint, Constraint)
        ]

        if dnachisel_constraints != []:
            if not DNACHISEL_AVAILABLE:
                raise ImportError("Spotted DNA Chisel constraints, while "
                                  "DNA Chisel is not installed.")
            canvas = DnaCanvas(sequence, dnachisel_constraints)
            constraints = [
                constraint for constraint in constraints
                if not isinstance(constraint, Constraint)
            ] + [lambda seq: canvas.all_constraints_pass()]

        return all(constraint(sequence) for constraint in constraints)


    def compute_supply_graph(self):
        """Return elements to plot the supply graph underlying this DnaSource

        Returns
        -------

        edges
          A list [(s1,s2), (s1,s3), (s2, s5)...] of couples of DnaSources in a
          supplier-supplied relationship.

        levels
          A list of lists [[s1,s2], [s4,s8,s9]...] of sources. The first
          sublist (first level) are all sources at the farthest distance from
          the current source in the supply graph, and the last sublist contains
          only the current DnaSource.
        """

        seen_sources = set()
        levels = defaultdict(lambda: [])
        edges = []

        def rec(source, depth=0):

            if source in seen_sources:
                return
            seen_sources.add(source)

            levels[depth].append(source)

            if hasattr(source, "dna_source"):
                edges.append((source.dna_source, source))
                rec(source.dna_source, depth + 1)

            if hasattr(source, "dna_sources"):
                for other in source.dna_sources:
                    edges.append((other, source))
                    rec(other, depth + 1)

        rec(self)
        levels = [levels[i] for i in sorted(levels.keys())][::-1]
        return edges, levels


class DnaSourcesComparator(DnaSource):
    """Special source that compares quotes from other DnaSources.

    Upon receiving a sequence, that source will submit the sequence to
    Downstream sources, which deliver each one optimal quote. The comparator
    then returns the one quote with the lowest price.

    Parameters
    ----------

    dna_sources
      List of `DnaSources` that

    memoize
      Whether the quotes should be kept in memory to avoid re-computing the
      same quote several times. Can accelerate computations but is
      RAM-expensive.

    """

    def __init__(self, dna_sources, memoize=False):
        self.dna_sources = dna_sources
        self.sequence_constraints = ()
        self.memoize = memoize
        self.memoize_dict = {}

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_ordering_plan
          If True, the ordering plan is added to the quote
        """
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
        self.extra_cost = assembly_method.cost
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
            location_filters=[
                lambda location: fl(sequence, location)
                for fl in assembly.location_filters
            ],
            segment_filters=[
                lambda (start, end): fl(sequence, start, end)
                for fl in assembly.segment_filters
            ],
            min_segment_length=assembly.min_segment_length,
            max_segment_length=assembly.max_segment_length,
            forced_cuts=assembly.force_cuts(sequence),
            initial_resolution=nucleotide_resolution,
            refine_resolution=refine_resolution,
            progress_bars=progress_bars,
            a_star_factor=a_star_factor,
            max_fragments=assembly.max_fragments
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
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_ordering_plan
          If True, the ordering plan is added to the quote
        """

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
        total_price = ordering_plan.total_price() + self.extra_cost
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

        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_ordering_plan
          If True, the ordering plan is added to the quote
        """

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
    """Class to represent databases of constructs which can be (in part) reused

    A blast database contains the sequences of all available constructs.
    Given a sequence, the PcrOutStation finds whether it is possible to order
    two primers to extract this sequence from the constructs in the BLAST
    database.

    Parameters
    ----------

    name
      Name of the PCR station (e.g. "Lab constructs PCR station")

    primers_dna_source
      DnaSource providing the primers (will typically be an ExternalDnaOffer)

    blast_database

    sequences

    pcr_homology_length

    max_overhang_length

    extra_cost

    extra_time

    max_amplicon_length

    blast_word_size

    memoize

    sequence_constraints

    """

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
                                     max_amplicon_length] +
                                    list(sequence_constraints))
        self.sequence_constraints = sequence_constraints
        self.memoize = memoize
        self.memoize_dict = {}
        self.sequences = sequences

    def get_hits(self, sequence):
        """Return the hits of the given sequence against the blast database
        (in Biopython format)"""
        if self.sequences is not None:
            result = []
            for dna_name, seq in self.sequences:
                match_coords = largest_common_substring(sequence, seq,
                                                        self.max_overhang_length)
                if match_coords:
                    result.append((dna_name, match_coords))
            return result
        else:
            record = blast_sequence(sequence, self.blast_database,
                                    perc_identity=100,
                                    word_size=self.blast_word_size)
            return [
                (alignment.hit_id, (hit.query_start, hit.query_end))
                for alignment in record.alignments
                for hit in alignment.hsps
            ]

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):
        """Return a price-optimal DnaQuote for the given sequence.

        It will find a possible hit in the blast database, find the primers to
        order for the PCR, compute the overall price and lead time, and return
        a quote.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_ordering_plan
          If True, the ordering plan is added to the quote
        """

        for subject, (hit_start, hit_end) in self.get_hits(sequence):

            largest_overhang = max(hit_start, len(sequence) - hit_end)

            if largest_overhang <= self.max_overhang_length:
                primer_l_end = hit_start + self.pcr_homology_length
                primer_left = sequence[:primer_l_end]
                primer_r_end = hit_end - self.pcr_homology_length
                primer_right = reverse_complement(sequence[primer_r_end:])

                primer_max_lead_time = (None if max_lead_time is None else
                                        max_lead_time - self.extra_time)
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
                        (0, primer_l_end): quotes[0],
                        (primer_r_end, len(sequence)): quotes[0]
                    })
                else:
                    ordering_plan = None
                return DnaQuote(self, sequence, accepted=True,
                                lead_time=overall_lead_time,
                                price=total_price,
                                ordering_plan=ordering_plan,
                                metadata={"subject": subject,
                                          "location": (hit_start, hit_end)})

        return DnaQuote(self, sequence, accepted=False,
                        message="No valid match found")

    def pre_blast(self, sequence):
        """Pre-compute the BLAST of the current sequence against the database.

        Once a pre-blast has been performed, this PcrOutStation becomes
        specialized on that sequence and its subsequences, do not feed it with
        another different sequence. Do `self.sequences=None` to reinitialize
        and de-specialize this PcrOutStation.

        Examples
        --------

        >>> pcr_station = PcrOutStation("some_blast_database")
        >>> top_station = # some assembly station depending on pcr_station
        >>> pcr_station.pre-blast(my_sequence)
        >>> top_station.get_quote(my_sequence)
        >>> pcr_station.sequences=None # de-specializes the pcr station.
        """
        self.sequences = None  # destroy current pre-blast (used by get_hits)
        self.sequences = [
            (subject, sequence[start:end])
            for subject, (start, end) in self.get_hits(sequence)
        ]


class PartsLibrary(DnaSource):
    """Class for collections of ready-to-assemble parts"""

    def __init__(self, name, parts_dict, memoize=False):
        self.name = name
        self.sequence_constraints = []
        self.parts_dict = parts_dict
        self.memoize = False

    def get_best_price(self, sequence, max_lead_time=None,
                       with_ordering_plan=False):
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_ordering_plan
          If True, the ordering plan is added to the quote
       """
        if sequence in self.parts_dict:
            return DnaQuote(self, sequence, accepted=True,
                            price=0, lead_time=0,
                            message="Part name: " + self.parts_dict[sequence])

        return DnaQuote(self, sequence, accepted=False,
                        message="Sequence not in the library")

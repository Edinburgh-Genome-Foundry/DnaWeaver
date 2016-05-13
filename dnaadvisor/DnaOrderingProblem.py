"""Definition and solution DNA ordering problems."""

from DnaOffer import DnaOfferEvaluation
from optimization import optimize_cuts_with_graph_twostep
from DnaOffer import DnaOrderingPlan
from copy import copy



class DnaOrderingProblem:
    """Class for the definition of ordering problems.

    In an ordering problem, we define a series of offers (a "pricing"
    proposed by a company and subject to constraints on the sequence)
    a DNA sequence, and an assembly method. ordering_problem.solve()


    Examples
    --------

    >>> problem = DnaOrderingProblem(sequence, offers=[company_1, company_2],
    >>>                              assembly_method=GibsonAssemblyMethod(20))
    >>> solution = problem.solve(min_segment_length=100,
    >>>                          max_segment_length=4000)
    >>> print (problem.ordering_plan_summary(solution))


    Parameters
    ----------

    sequence
      A "ATGC" string representing a DNA sequence

    offers
      A list of DnaOffer objects

    assembly_method
      An AssemblyMethod object - assembly method that will be used.

    location_filters
      A list or tuple of functions f(index)->boolean used for prefiltering
      the possible cut locations. Only the locations verifying f(index)==True
      for all f in location_filters will be considered in the problem.

    segment_filters
       A list or tuple of functions f((start, end))->boolean used for
       prefiltering the segments of the sequence. Only the segments verifying
       f((start, end))==True for all f in segments_filters will be considered
       in the problem.

     cuts_number_penalty
       A number representing a penalty on the numbers of cuts. A large number
       leads to solutions with fewer fragments (but more expensive than the
       optimal solution). The score optimized is total_price + P*Nfragments
       where P is the cuts_number_penalty.

    """

    def __init__(self, sequence, offers, assembly_method,
                 location_filters=(), segment_filters=(),
                 cuts_number_penalty=0):
        self.sequence = sequence
        self.offers = offers
        self.assembly_method = assembly_method
        self.segment_filters = segment_filters
        self.location_filters = location_filters
        self.cuts_number_penalty = cuts_number_penalty

    def offers_for_segment(self, segment):
        fragment = self.assembly_method.compute_fragment_sequence(
            segment, self.sequence
        )
        result = []
        for offer in self.offers:
            evaluation = DnaOfferEvaluation(fragment, offer, segment)
            if evaluation.is_orderable:
                result.append(evaluation)
        return result

    def best_price_for_segment(self, segment):
        """Return the best price for a segment.

        First a sequence is computed by extracting the segment from the
        problem's sequence, and possibly adding flanking overhangs (depending
        on the assembly method used, e.g. Gibson Assembly).

        Then the minimal price is computed among all the problem's DNA offers
        whose constraints are verified by the sequence.
        """
        offers = self.offers_for_segment(segment)
        if offers == []:
            return -1
        else:
            return min(offer.price for offer in offers)

    def compute_fragments_sequences(self, cuts):
        """Return all fragments necessary for assembly, depending on the cuts.
        """
        return self.assembly_method.compute_fragments_sequences(
            cuts, self.sequence
        )

    def compute_all_offers(self, cuts):
        """Return all valid offers for all fragments corresponding to the cuts.
        """
        segments_fragments = self.compute_fragments_sequences(cuts)
        offers = {fragment: [] for fragment in segments_fragments.values()}
        for segment, fragment in segments_fragments.items():
            for offer in self.offers:
                evaluation = DnaOfferEvaluation(fragment, offer, segment)
                if evaluation.is_orderable:
                    offers[fragment].append(evaluation)
        return offers

    def find_best_ordering_plan(self, cuts):
        """Find the best price-wise ordering plan corresponding to the cuts.

        First the DNA fragments corresponding to the cuts are computed.
        For each fragment, we evaluate all offer and keep the valid offer with
        the lowest price.

        Parameters
        ----------

        cuts
          A list of integers between 0 and len(self.sequence) representing
          cut locations.

        Returns
        -------

        ordering_plan
          A DnaOrderingPlan object

        """
        offers = self.compute_all_offers(cuts)
        best_offers = {}
        for fragment, fragment_offers in offers.items():
            if fragment_offers == []:
                best_offers[fragment] = None
            else:
                best_offer = min(fragment_offers, key=lambda o: o.price)
                best_offers[fragment] = best_offer
        return DnaOrderingPlan(best_offers.values(),
                               full_sequence=self.sequence)

    def solve(self, min_segment_length=500, max_segment_length=2000,
              nucleotide_resolution=1, refine_resolution=1,
              forced_cuts=None):
        """Solve the ordering problem.

        Examples
        --------

        >>> solution = ordering_problem.solve()
        >>> print (ordering_problem.ordering_plan_summary(solution))

        Parameters
        ----------

        min_segment_length (int)
          Minimal length of the segments (without possible overhangs)

        max_segment_length (int)
          Maximal length of the segments (without possible overhangs)

        nucleotide_resolution
          Resolution to use in the solution of the problem (before refinement).
          If equal to 10, only every 10th nucleotide will be considered as a
          possible cutting point (thus reducing the size of the problem)

        refine_resolution
          Resolution to apply during the refinement.

        Returns
        -------

        ordering_plan
          A DnaOrderingPlan object summarizing what sequences to order, and
          from which offers they should be ordered.

        """

        if forced_cuts is not None:
            forced_cuts = list(sorted(forced_cuts))
            fragments =self.compute_fragments_sequences(forced_cuts)
            offers = []
            for segment in sorted(fragments.keys()):
                sub_problem = copy(self)
                sub_problem.sequence = fragments[segment]
                sub_offers = sub_problem.solve(
                    min_segment_length=min_segment_length,
                    max_segment_length=max_segment_length,
                    nucleotide_resolution=nucleotide_resolution,
                    refine_resolution=refine_resolution
                ).plan
                for offer in sub_offers:
                    start, stop = offer.segment
                    offer.segment = (start + segment[0],
                                     min(segment[1], stop + segment[0]))
                    offers.append(offer)
            return DnaOrderingPlan(offers, full_sequence=self.sequence)

        best_cuts = optimize_cuts_with_graph_twostep(
            sequence_length=len(self.sequence),
            segment_score_function=self.best_price_for_segment,
            cuts_number_penalty=self.cuts_number_penalty,
            location_filters=self.location_filters,
            segment_filters=self.segment_filters,
            initial_resolution=nucleotide_resolution,
            min_segment_length=min_segment_length,
            max_segment_length=max_segment_length,
            refine_resolution=refine_resolution
        )

        return self.find_best_ordering_plan(best_cuts)

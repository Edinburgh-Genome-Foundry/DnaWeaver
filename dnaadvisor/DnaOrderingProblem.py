"""Definition of a problem."""

from DnaOffer import DnaOfferEvaluation
from optimization import optimize_cuts_with_graph_twostep


class DnaOrderingProblem:

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
        offers = self.offers_for_segment(segment)
        if offers == []:
            return -1
        else:
            return min(offer.price for offer in offers)

    def compute_all_offers(self, cuts):
        segments_fragments = self.compute_fragments_sequences(cuts)
        offers = {fragment: [] for fragment in segments_fragments.values()}
        for segment, fragment in segments_fragments.items():
            for offer in self.offers:
                evaluation = DnaOfferEvaluation(fragment, offer, segment)
                if evaluation.is_orderable:
                    offers[fragment].append(evaluation)
        return offers

    def compute_fragments_sequences(self, cuts):
        return self.assembly_method.compute_fragments_sequences(
            cuts, self.sequence
        )

    def find_best_ordering_plan(self, cuts):
        offers = self.compute_all_offers(cuts)
        best_offers = {}
        for fragment, fragment_offers in offers.items():
            if fragment_offers == []:
                best_offers[fragment] = None
            else:
                best_offer = min(fragment_offers, key=lambda o: o.price)
                best_offers[fragment] = best_offer
        return best_offers

    @staticmethod
    def ordering_plan_summary(plan):
        evaluations = sorted(plan.values(), key=lambda o: o.segment)
        plan = "\n  ".join(str(e) for e in evaluations)
        total_price = sum(ev.price for ev in evaluations)
        return "Ordering plan:\n  %s\n  Total price: %d" % (plan, total_price)

    def solve(self, min_segment_length=500, max_segment_length=2000,
              nucleotide_resolution=1, refine_resolution=1):

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

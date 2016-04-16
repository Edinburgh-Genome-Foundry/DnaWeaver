"""Definition of a problem."""
from dnachisel import DnaCanvas

class OfferEvaluation:

    def __init__(self, sequence, offer, zone=None):

        self.sequence = sequence
        self.offer = offer
        self.zone = zone
        canvas = DnaCanvas(sequence)
        if not all([c.evaluate(canvas).passes for c in offer.constraints]):
            self.is_orderable = False
            self.price = None
        else:
            self.is_orderable = True
            self.price = offer.pricing(sequence)

    def __repr__(self):
        if not self.is_orderable:
            return "%s (not orderable)" % self.offer.name
        else:
            return "%s%s %.02f$" % ("" if self.zone is None else ("%s " % str(self.zone)),
                                    self.offer.name, self.price)


class DnaOffer:

    def __init__(self, name, constraints, pricing):
        self.name = name
        self.constraints = constraints
        self.pricing = pricing

    def __repr__(self):
        return self.name


class DnaOrderingProblem:

    def __init__(self, sequence, offers, assembly_method):
        self.sequence = sequence
        self.offers = offers
        self.assembly_method = assembly_method

    def compute_fragments_zones(self, cuts):
        L = len(self.sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        return zip(cuts, cuts[1:])

    def validate_cuts(self, cuts):
        return self.assembly_method.validate_cuts(cuts, self.sequence)

    def compute_fragments_sequences(self, cuts):
        return self.assembly_method.compute_fragments_sequences(
            cuts, self.sequence
        )

    def compute_all_offers(self, cuts):
        zones_fragments = self.compute_fragments_sequences(cuts)
        offers = {fragment: [] for fragment in zones_fragments.values()}
        for zone, fragment in zones_fragments.items():
            for offer in self.offers:
                evaluation = OfferEvaluation(fragment, offer, zone)
                if evaluation.is_orderable:
                    offers[fragment].append(evaluation)
        return offers

    def compute_orderability_and_costs(self, cuts):
        fragments = self.compute_fragments_sequences(cuts)
        offers = {fragment: [] for fragment in fragments}
        for fragment in fragments:
            for offer in self.offers:
                evaluation = offer.evaluate(fragment)
                if evaluation.is_orderable:
                    offers[fragment].append(evaluation)

    def score_cuts(self, cuts):
        INVALID_CUTS = 1000000.0
        if not self.validate_cuts(cuts):
            return INVALID_CUTS

        offers = self.compute_all_offers(cuts)
        if any(offers_list == [] for offers_list in offers.values()):
            return INVALID_CUTS
        best_price = sum(min([ev.price for ev in offers_list])
                         for offers_list in offers.values())
        return best_price

    def find_best_offers(self, cuts):
        offers = self.compute_all_offers(cuts)
        best_offers = {}
        for fragment, fragment_offers in offers.items():
            if fragment_offers == []:
                best_offers[fragment] = None
            else:
                best_offer = min(fragment_offers, key=lambda o: o.price)
                best_offers[fragment] = best_offer
        return best_offers


class AssemblyMethod:
    pass


class GibsonAssemblyMethod(AssemblyMethod):

    def __init__(self, overlap=40):
        self.overlap = overlap

    def compute_fragments_sequences(self, cuts, sequence):
        L = len(sequence)
        cuts = sorted(list(set([0, L] + cuts)))
        zones = zip(cuts, cuts[1:])
        overlap = self.overlap
        return {
            zone: sequence[max(0, zone[0] - overlap):
                           min(L, zone[1] + overlap)]
            for zone in zones
        }

    def validate_cuts(self, cuts, sequence):
        return True

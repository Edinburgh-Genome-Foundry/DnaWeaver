from dnachisel import DnaCanvas


class DnaOfferEvaluation:

    def __init__(self, sequence, offer, segment=None):

        self.sequence = sequence
        self.offer = offer
        self.segment = segment
        canvas = DnaCanvas(sequence, constraints=offer.constraints)
        if not canvas.all_constraints_pass():
            self.is_orderable = False
            self.price = None
        else:
            self.is_orderable = True
            self.price = offer.pricing(sequence)

    def __repr__(self):
        if not self.is_orderable:
            return "%s (not orderable)" % self.offer.name
        else:
            return "%s%s %.02f$" % ("" if self.segment is None
                                    else ("%s " % str(self.segment)),
                                    self.offer.name, self.price)


class DnaOffer:

    def __init__(self, name, constraints, pricing):
        self.name = name
        self.constraints = constraints
        self.pricing = pricing

    def __repr__(self):
        return self.name

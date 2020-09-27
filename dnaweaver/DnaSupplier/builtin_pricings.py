"""Pricings for suppliers such as CommercialDnaOffer."""


class PerBasepairPricing:
    """Class for pricings of the form X/basepair + fixed cost.

    Parameters
    ----------

    per_basepair_pricing
      Variable cost per basepair on top of the fixed cost.

    fixed_cost
      Fixed cost added to the variable cost for each sequence.

    Examples
    --------

    >>> commercial_offer = CommercialDnaOffer(
    >>>     name="Company 1",
    >>>     pricing=PerBasepairPricing(0.14, fixed_cost=30)
    >>> )
    """

    def __init__(self, per_basepair_price, fixed_cost=0):
        self.per_basepair_price = per_basepair_price
        self.min_basepair_price = per_basepair_price
        self.fixed_cost = fixed_cost

    def __call__(self, sequence):
        return len(sequence) * self.per_basepair_price + self.fixed_cost

    def __str__(self):
        if self.fixed_cost != 0:
            return "$(%.02f + %.03f/bp)" % (self.fixed_cost, self.per_basepair_price,)
        return "$%.03f/bp" % self.per_basepair_price


class FixedCostPricing:
    """Class for pricings that return the same price for each sequence.

    Examples
    --------

    >>> commercial_offer = CommercialDnaOffer(
    >>>     name="Company 1",
    >>>     pricing=FixedCostPricing(5)
    >>> )
    """

    def __init__(self, fixed_cost):
        self.fixed_cost = fixed_cost

    def __call__(self, sequence):
        return self.fixed_cost

    def __str__(self):
        return "$%.03f/order" % self.fixed_cost

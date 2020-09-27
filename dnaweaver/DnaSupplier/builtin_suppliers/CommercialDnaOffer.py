import numpy as np
from ...tools import functions_list_to_string
from ...DnaQuote import DnaQuote
from ..DnaSupplier import DnaSupplier


from Bio.Restriction import AllEnzymes

from ..builtin_constraints import (
    GcContentConstraint,
    NoPatternConstraint,
    SequenceLengthConstraint,
)
from ..builtin_pricings import PerBasepairPricing


class CommercialDnaOffer(DnaSupplier):
    """External/Commercial source of DNA."""

    class_description = "External DNA offer"
    operation_type = "order"
    report_fa_symbol = u""
    report_fa_symbol_plain = "shopping-cart"
    report_color = "#ffeeee"

    def __init__(
        self, name, pricing, sequence_constraints=(), lead_time=0, memoize=False,
    ):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.pricing = pricing
        self.lead_time = lead_time if callable(lead_time) else (lambda *a: lead_time)
        self.memoize = memoize
        self.memoize_dict = {}
        self.min_basepair_price = 0
        if hasattr(pricing, "min_basepair_price"):
            self.min_basepair_price = pricing.min_basepair_price

        if hasattr(pricing, "fixed_cost"):
            max_length = min(
                [
                    cst.max_length
                    for cst in sequence_constraints
                    if isinstance(cst, SequenceLengthConstraint)
                ],
                default=None,
            )
            if max_length is not None:
                self.min_basepair_price += pricing.fixed_cost / max_length

    def __repr__(self):
        return self.name

    def get_best_price(self, sequence, max_lead_time=None, with_assembly_plan=False):
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quote.

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_assembly_plan
          If True, the assembly plan is added to the quote.
        """

        lead_time = self.lead_time(sequence)
        price = self.pricing(sequence)
        accepted = (max_lead_time == np.inf) or (lead_time <= max_lead_time)
        return DnaQuote(
            self, sequence, accepted=accepted, lead_time=lead_time, price=price
        )

    def get_best_lead_time_under_price_limit(self, sequence, max_price):

        lead_time = self.lead_time(sequence)
        price = self.pricing(sequence)
        return DnaQuote(
            self,
            sequence,
            price=price,
            lead_time=lead_time,
            accepted=price < max_price,
        )

    def additional_dict_description(self):
        constraints = functions_list_to_string(self.sequence_constraints)
        return {
            "price function": functions_list_to_string([self.pricing]),
            "sequence constraints": constraints,
        }

    @staticmethod
    def from_dict(data):
        constraints = []
        if "gc_range" in data:
            mini, maxi = data["gc_range"]
            if (mini, maxi) != (0, 100):
                constraints.append(GcContentConstraint(0.01 * mini, 0.01 * maxi))
        if "size_range" in data:
            mini, maxi = data["size_range"]
            constraints.append(SequenceLengthConstraint(mini, maxi))
        if "forbidden" in data:
            ENZYMES = list(str(e) for e in AllEnzymes)
            for pattern in data["forbidden"].split(", "):
                if pattern == "":
                    continue
                enzyme = None
                if pattern in ENZYMES:
                    enzyme = pattern
                    pattern = None
                constraints.append(NoPatternConstraint(pattern=pattern, enzyme=enzyme))

        return CommercialDnaOffer(
            name=data["name"],
            pricing=PerBasepairPricing(data["cost_per_bp"]),
            lead_time=data["lead_time"],
            memoize=data.get("memoize", False),
            sequence_constraints=constraints,
        )

from ..DnaQuote import DnaQuote
from .DnaSource import DnaSource
from ..tools import functions_list_to_string

import numpy as np

class CommercialDnaOffer(DnaSource):
    """External/Commercial source of DNA"""

    class_description = "External DNA offer"
    operation_type = "order"
    report_fa_symbol = u""
    report_fa_symbol_plain = "shopping-cart"
    report_color = "#ffeeee"

    def __init__(self, name, pricing, sequence_constraints=(),
                 lead_time=None, memoize=False):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.pricing = pricing
        self.lead_time = (lead_time if callable(lead_time)
                          else (lambda *a: lead_time))
        self.memoize = memoize
        self.memoize_dict = {}
        if hasattr(pricing, 'min_basepair_price'):
            self.min_basepair_price = pricing.min_basepair_price

    def __repr__(self):
        return self.name

    def get_best_price(self, sequence, max_lead_time=None,
                       with_assembly_plan=False):
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_assembly_plan
          If True, the assembly plan is added to the quote
        """

        lead_time = self.lead_time(sequence)
        price = self.pricing(sequence)
        accepted = (max_lead_time == np.inf) or (lead_time <= max_lead_time)
        return DnaQuote(self, sequence, accepted=accepted,
                        lead_time=lead_time, price=price)

    def get_best_lead_time_under_price_limit(self, sequence, max_price):

        lead_time = self.lead_time(sequence)
        price = self.pricing(sequence)
        return DnaQuote(self, sequence, price=price, lead_time=lead_time,
                        accepted=price < max_price)

    def additional_dict_description(self):
        constraints = functions_list_to_string(self.sequence_constraints)
        return {
            "price function": functions_list_to_string([self.pricing]),
            "sequence constraints": constraints
        }

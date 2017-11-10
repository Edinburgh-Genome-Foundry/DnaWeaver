from .DnaSource import DnaSource
from ..DnaQuote import DnaQuote

class PartsLibrary(DnaSource):
    """Class for collections of ready-to-assemble parts.


    This class is admitedly under-developed and could be expanded-subclassed
    to accomodate with the different kinds of registries etc.
    """
    class_description = "Parts Library"
    operation_type = "library"
    report_fa_symbol = u""
    report_fa_symbol_plain = "book"

    report_color = "#feeefe"

    def __init__(self, name, parts_dict, memoize=False,
                 sequence_constraints=()):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.parts_dict = parts_dict
        self.memoize = False

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
        if sequence in self.parts_dict:
            return DnaQuote(self, sequence, accepted=True,
                            price=0, lead_time=0,
                            message="Part name: " + self.parts_dict[sequence])

        return DnaQuote(self, sequence, accepted=False,
                        message="Sequence not in the library")

    def additional_dict_description(self):
        return {
            "flanks length": self.flanks_length
        }


class GoldenGatePartsLibrary(PartsLibrary):
    """Library of parts for Golden Gate Assembly"""
    class_description = "Golden Gate parts library"

    def __init__(self, name, parts_dict, flanks_length=7, memoize=False,
                 sequence_constraints=()):
        self.name = name
        self.sequence_constraints = sequence_constraints
        self.parts_dict = parts_dict
        self.flanks_length = flanks_length
        self.memoize = False

    def suggest_cuts(self, sequence):
        suggested_cuts = []
        # + 2 is because the cut falls in the middle of the 4bp linker:
        flank = self.flanks_length + 2
        for part in self.parts_dict:
            segment = part[flank:-flank]
            i = sequence.find(segment)
            if i != -1:
                suggested_cuts += [i, i + len(segment)]
        return sorted(list(set(suggested_cuts)))

    def additional_dict_description(self):
        return {
            "class": "Golden Gate parts library",
            "operation_type": "library",
            "flanks length": self.flanks_length
        }

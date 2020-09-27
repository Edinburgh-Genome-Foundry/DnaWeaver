from ..DnaSupplier import DnaSupplier
from ...DnaQuote import DnaQuote


class PcrLinearizationStation(DnaSupplier):
    """Station which will linearize another station's (circular) output via PCR.
    
    Parameters
    ----------
    """

    class_description = "Linearizes a sequence for downstream assembly"
    report_fa_symbol = u""
    report_fa_symbol_plain = "exchange"
    report_color = "#eefefe"
    operation_type = "PCR"

    def __init__(
        self,
        name,
        supplier,
        primers_supplier,
        homology_selector,
        sequence_constraints=(),
        extra_time=0,
        extra_cost=0,
        memoize=False,
    ):
        self.name = name
        self.supplier = supplier
        self.primers_supplier = primers_supplier
        # only for network reconstitution
        self.suppliers = [supplier, primers_supplier]
        self.homology_selector = homology_selector
        self.sequence_constraints = sequence_constraints
        self.extra_time = extra_time
        self.extra_cost = extra_cost
        self.memoize = memoize
        self.min_basepair_price = supplier.min_basepair_price

    def get_best_price(
        self, sequence, max_lead_time=None, with_assembly_plan=False,
    ):
        """Return a price-optimal DnaQuote for the given sequence.

        It will find a possible hit in the blast database, find the primers to
        order for the PCR, compute the overall price and lead time, and return
        a quote.

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
        if max_lead_time is not None:
            suppliers_max_lead_time = max_lead_time - self.extra_time
        else:
            suppliers_max_lead_time = None
        supplier_quote = self.supplier.get_quote(
            sequence, max_lead_time=suppliers_max_lead_time
        )
        if not supplier_quote.accepted:
            return DnaQuote(
                self, sequence, accepted=False, message="Refused by supplier"
            )

        left_location = self.homology_selector.compute_segment_location(
            sequence=sequence, index=0
        )
        right_location = self.homology_selector.compute_segment_location(
            sequence=sequence, index=len(sequence)
        )
        if (left_location is None) or (right_location is None):
            return DnaQuote(
                self,
                sequence,
                accepted=False,
                message="Could not design a suitable primer",
            )
        _, left_end = left_location
        right_start, _ = right_location

        left_primer = sequence[:left_end]
        right_primer = sequence[right_start:]

        primers_quotes = [
            self.primers_supplier.get_quote(
                primer, max_lead_time=suppliers_max_lead_time
            )
            for primer in [left_primer, right_primer]
        ]
        if not all(quote.accepted for quote in primers_quotes):
            return DnaQuote(
                self,
                sequence,
                accepted=False,
                message="Primers refused by primers supplier",
            )
        all_quotes = [supplier_quote] + primers_quotes
        if max_lead_time is not None:
            overall_lead_time = (
                max(quote.lead_time for quote in all_quotes) + self.extra_time
            )
        else:
            overall_lead_time = None
        total_price = sum(quote.price for quote in all_quotes) + self.extra_cost

        if with_assembly_plan:
            assembly_plan = {
                (0, len(sequence)): all_quotes[0],
                (0, int(left_end)): all_quotes[1],
                (int(right_start), len(sequence)): all_quotes[2],
            }
        else:
            assembly_plan = None

        return DnaQuote(
            self,
            sequence,
            accepted=True,
            lead_time=overall_lead_time,
            price=total_price,
            assembly_plan=assembly_plan,
            metadata={"location": (0, len(sequence)),},
        )

    def additional_dict_description(self):
        return {
            "primers DNA source": self.primers_supplier.name,
            "primers homology selector": str(self.homology_selector),
        }

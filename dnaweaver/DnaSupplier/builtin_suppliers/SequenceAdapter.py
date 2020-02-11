from ..DnaSupplier import DnaSupplier


class SequenceAdapter(DnaSupplier):
    """Transparent source adding sequence flanks between suppliers and clients.
    """

    class_description = "Adapts a sequence between a supplier and its clients"
    operation_type = "sequence extension"
    report_fa_symbol = u"\uf07e"
    report_fa_symbol_plain = "arrows-h"
    report_color = "#ffeeee"

    def __init__(
        self,
        name,
        supplier,
        left_addition="",
        right_addition="",
        memoize=False,
    ):
        self.name = name
        self.supplier = supplier
        self.suppliers = [supplier] # only for network reconstitution
        self.left_addition = left_addition
        self.right_addition = right_addition
        self.sequence_constraints = ()
        self.memoize = memoize
        self.min_basepair_price = supplier.min_basepair_price
        if hasattr(supplier, "min_basepair_price"):
            self.min_basepair_price = supplier.min_basepair_price

    def __repr__(self):
        return self.name

    def get_best_price(
        self,
        sequence,
        max_lead_time=None,
        with_assembly_plan=False,
    ):
        extended_sequence = self.left_addition + sequence + self.right_addition
        return self.supplier.get_quote(
            extended_sequence,
            max_lead_time=max_lead_time,
            with_assembly_plan=with_assembly_plan,
        )

    @staticmethod
    def from_dict(data):
        return SequenceAdapter(
            name=data["name"],
            supplier=data["supplier"],
            left_addition=data["left_addition"],
            right_addition=data["right_addition"],
        )

    def suggest_cuts(self, sequence):
        return self.supplier.suggest_cuts(sequence)

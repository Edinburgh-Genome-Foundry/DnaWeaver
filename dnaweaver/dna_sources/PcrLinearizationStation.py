from .DnaSource import DnaSource


class PcrLinearizationStation(DnaSource):
    """PCR-Out a fragment from a vector to linearize it for use in subsequent
    assemblies such as Gibson assembly."""

    report_fa_symbol = u""
    report_fa_symbol_plain = "exchange"
    report_color = "#eefefe"
    operation_type = "PCR"

    def __init__(
        self,
        supplier,
        primers_supplier,
        primer_melting_temperature=50,
        sequence_constraints=(),
    ):
        self.supplier = supplier
        self.primers_supplier = primers_supplier
        self.primer_melting_temperature = primer_melting_temperature
        self.sequence_constraints = sequence_constraints

    def additional_dict_description(self):
        return {
            "primers DNA source": self.primers_supplier.name,
            "primers melting temp.": self.primer_melting_temperature,
        }

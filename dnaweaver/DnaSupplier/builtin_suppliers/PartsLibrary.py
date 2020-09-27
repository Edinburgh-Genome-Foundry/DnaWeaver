from Bio import SeqIO
from ...DnaQuote import DnaQuote
from ..DnaSupplier import DnaSupplier


class PartsLibrary(DnaSupplier):
    """Class for collections of ready-to-assemble parts.

    This class is admittedly under-developed and could be expanded-subclassed
    to accommodate the different kinds of registries etc.
    """

    class_description = "Parts Library"
    operation_type = "library"
    report_fa_symbol = u""
    report_fa_symbol_plain = "book"

    report_color = "#feeefe"
    collections_by_id = {}
    library_classes = {}

    def __init__(
        self,
        name,
        parts_dict=None,
        fasta_file=None,
        memoize=False,
        price_per_part=0,
        lead_time=0,
        sequence_constraints=(),
    ):
        self.name = name
        self.price_per_part = price_per_part
        self.lead_time = lead_time
        self.sequence_constraints = sequence_constraints
        if fasta_file is not None:
            parts_dict = {
                record.id: str(record.seq).upper()
                for record in SeqIO.parse(fasta_file, "fasta")
            }
        self.parts_dict = parts_dict
        self.inverted_parts_dict = {v: k for k, v in parts_dict.items()}
        self.sequences_set = set(self.inverted_parts_dict)
        self.memoize = memoize
        self.memoize_dict = {}

    def get_best_price(
        self, sequence, max_lead_time=None, with_assembly_plan=False,
    ):
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
        sequence = self.preprocess_sequence(sequence)
        if sequence in self.sequences_set:
            part_name = self.inverted_parts_dict[sequence]
            return DnaQuote(
                self,
                sequence,
                accepted=True,
                price=self.price_per_part,
                lead_time=self.lead_time,
                message="Part: " + part_name,
                metadata={"part_name": part_name},
            )

        return DnaQuote(
            self, sequence, accepted=False, message="Sequence not in the library",
        )

    def preprocess_sequence(self, sequence):
        """Can be used by subclasses e.g. to anonymize wildcard nucleotides"""
        return sequence

    def additional_dict_description(self):
        return {"flanks length": self.flanks_length}

    @classmethod
    def from_dict(cls, data):
        parameters = cls.collections_by_id[data["collection"]]
        library_class = parameters.pop("library_class")
        if library_class in cls.library_classes:
            library_class = cls.library_classes[library_class]

        def get(param, default):
            return data.get(param, parameters.get(param, default))

        return library_class(
            name=get("name", "library"),
            price_per_part=get("price_per_part", 0),
            lead_time=get("lead_time", 0),
            parts_dict=get("parts_dict", None),
            fasta_file=get("fasta_file", None),
            memoize=get("memoize", None),
        )


class GoldenGatePartsLibrary(PartsLibrary):
    """Library of parts for Golden Gate Assembly."""

    class_description = "Golden Gate parts library"

    def __init__(
        self,
        name,
        parts_dict=None,
        fasta_file=None,
        price_per_part=0,
        lead_time=0,
        flanks_length=7,
        memoize=False,
        sequence_constraints=(),
    ):
        PartsLibrary.__init__(
            self,
            name,
            parts_dict=parts_dict,
            fasta_file=fasta_file,
            memoize=memoize,
            sequence_constraints=sequence_constraints,
        )
        self.flanks_length = flanks_length

    def suggest_cuts(self, sequence):
        suggested_cuts = []
        # + 2 is because the cut falls in the middle of the 4bp linker:
        flank = self.flanks_length
        for part, part_sequence in self.parts_dict.items():
            segment = part_sequence[flank:-flank]
            i = sequence.find(segment)
            if i != -1:
                suggested_cuts += [i + 2, i + len(segment) - 2]
        return sorted(list(set(suggested_cuts)))

    def suggest_segments(self, sequence):
        suggested_segments = []
        # + 2 is because the cut falls in the middle of the 4bp linker:
        flank = self.flanks_length
        for part, part_sequence in self.parts_dict.items():
            segment = part_sequence[flank:-flank]
            i = sequence.find(segment)
            if i != -1:
                L = len(segment)
                suggested_segments.append(((i + 2, i + L - 2), part))
        return sorted(set(suggested_segments))

    @classmethod
    def preprocess_sequence(cls, sequence):
        """Can be used by subclasses e.g. to anonymize wildcard nucleotides"""
        return sequence[:6] + "N" + sequence[7:-7] + "N" + sequence[-6:]

    def additional_dict_description(self):
        return {
            "class": "Golden Gate parts library",
            "operation_type": "library",
            "flanks length": self.flanks_length,
        }


PartsLibrary.library_classes.update(
    {"library": PartsLibrary, "golden_gate": GoldenGatePartsLibrary}
)

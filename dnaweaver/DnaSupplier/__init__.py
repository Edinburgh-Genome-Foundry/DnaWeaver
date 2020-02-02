"""DnaWeaver Reports implements reporting methods for Dna Weaver outputs."""

from .DnaSupplier import DnaSupplier
from .builtin_suppliers import (
    SequenceAdapter,
    CommercialDnaOffer,
    DnaAssemblyStation,
    DnaSuppliersComparator,
    PartsLibrary,
    GoldenGatePartsLibrary,
    PcrExtractionStation,
    PcrLinearizationStation,
)

DnaSupplier.default_suppliers_dict = {
    "commercial": CommercialDnaOffer,
    "assembly": DnaAssemblyStation,
    "library": PartsLibrary,
    "golden_gate_library": GoldenGatePartsLibrary,
    "pcr": PcrExtractionStation,
    "comparator": DnaSuppliersComparator,
    "main": DnaSuppliersComparator,
}

__all__ = [
    "DnaSupplier",
    "SequenceAdapter",
    "CommercialDnaOffer",
    "DnaAssemblyStation",
    "DnaSuppliersComparator",
    "PartsLibrary",
    "GoldenGatePartsLibrary",
    "PcrExtractionStation",
    "PcrLinearizationStation",
]

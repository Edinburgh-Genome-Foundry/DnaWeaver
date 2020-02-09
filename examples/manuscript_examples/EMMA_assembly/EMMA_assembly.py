import os
from Bio import SeqIO
from dnaweaver import (
    CommercialDnaOffer,
    GoldenGateAssemblyMethod,
    FixedCostPricing,
    PerBasepairPricing,
    DnaAssemblyStation,
    GoldenGatePartsLibrary,
    SequenceLengthConstraint,
    PcrExtractionStation,
    TmSegmentSelector,
    load_record,
)


# DEFINE THE SUPPLY NETWORK

# The EMMA collection of mammalian genetic parts
emma_collection = GoldenGatePartsLibrary(
    "EMMA",
    parts_dict={
        record.id: str(record.seq)
        for record in SeqIO.parse("emma_parts.fa", "fasta")
    },
    memoize=True,
)

# A medium-price vendor who can provide long parts
company_ingen = CommercialDnaOffer(
    name="DeluxeDNA",
    pricing=PerBasepairPricing(0.14),
    sequence_constraints=[SequenceLengthConstraint(max_length=2400)],
)

# A cheap vendor who can provide small parts < 1kb
company_tdi = CommercialDnaOffer(
    name="CheapDNA",
    pricing=PerBasepairPricing(0.08),
    sequence_constraints=[SequenceLengthConstraint(max_length=1000)],
)

# An oligos vendor (for oligo assembly and )
company_oligo = CommercialDnaOffer(
    name="Oligo vendor",
    pricing=FixedCostPricing(5),
    sequence_constraints=[SequenceLengthConstraint(max_length=100)],
)


mouse_pcr_station = PcrExtractionStation(
    name="E. coli",
    extra_cost=10,
    homology_selector=TmSegmentSelector(),
    primers_supplier=company_oligo,
    blast_database=os.path.join("..", "..", "data", "ecoli_blast_db", "ecoli"),
    memoize=True,
)

assembly_station = DnaAssemblyStation(
    name="Golden Gate Station",
    assembly_method=GoldenGateAssemblyMethod(
        min_segment_length=40, max_segment_length=5000, enzyme="BsmBI"
    ),
    supplier=[
        company_ingen,
        emma_collection,
        mouse_pcr_station,
        company_tdi,
    ],
    coarse_grain=100,
    fine_grain=None,
    logger="bar",
)

for station, color in [
    (company_ingen, "#fdfdb3"),
    (company_oligo, "#b8b6e6"),
    (company_tdi, "#c8f8c4"),
    (mouse_pcr_station, "#eff5f5"),
    (emma_collection, "#f8c3c3"),
]:
    station.report_color = color


# LOAD THE SEQUENCE AND PREPARE THE NETWORK FOR IT. HERE THE
# PREPARATION WILL ONLY MAKE THE PCR STATION PRE-BLAST THE SEQUENCE

record = load_record("example_sequence.gb")
assembly_station.prepare_network_on_sequence(record)
quote = assembly_station.get_quote(record, with_assembly_plan=True)
assembly_report = quote.to_assembly_plan_report()
assembly_report.write_full_report("EMMA_assembly_report")
print("Success! See figure assembly_plans.png.")


import os
from dnaweaver import (CommercialDnaOffer, GoldenGateAssemblyMethod,
                       PerBasepairPricing, DnaAssemblyStation,
                       GoldenGatePartsLibrary, SequenceLengthConstraint,
                       load_record)
THIS_DIR = os.path.join("tests", "test_scenarios",
                        "test_emma_homology_cassette")
EMMA_PATH = os.path.join(THIS_DIR, "emma_parts.fa")
SEQUENCE_PATH = os.path.join(THIS_DIR, "emma_sequence.gb")

def test_emma_construct():
    emma_collection = GoldenGatePartsLibrary("EMMA", fasta_file=EMMA_PATH,
                                             memoize=True)
    company_ingen = CommercialDnaOffer(
        name="InGen",
        pricing=PerBasepairPricing(0.14),
        sequence_constraints=[SequenceLengthConstraint(max_length=2400)]
    )

    company_tdi = CommercialDnaOffer(
        name="TDI",
        pricing=PerBasepairPricing(0.08),
        sequence_constraints=[SequenceLengthConstraint(max_length=600)]
    )

    assembly_station = DnaAssemblyStation(
        name='GoldenGate Assembly Station',
        assembly_method=GoldenGateAssemblyMethod(
            min_segment_length=40,
            max_segment_length=5000,
            enzyme='BsmBI',
    #         max_fragments=8
        ),
        supplier=[company_ingen, emma_collection, company_tdi],
        coarse_grain=100,
        fine_grain=10,
        logger='bar',
        a_star_factor='auto'
    )
    record = load_record(SEQUENCE_PATH)
    sequence = str(record.seq)
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
    emma_parts = [p for p in list(quote.assembly_plan.values())
              if p.source.name == 'EMMA']
    assert len(emma_parts) == 6

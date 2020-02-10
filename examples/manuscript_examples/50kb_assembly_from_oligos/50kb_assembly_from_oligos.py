from dnaweaver import (
    CommercialDnaOffer,
    DnaAssemblyStation,
    GibsonAssemblyMethod,
    OligoAssemblyMethod,
    TmSegmentSelector,
    FixedSizeSegmentSelector,
    PerBasepairPricing,
    SequenceLengthConstraint,
)

# OLIGO COMPANY

oligo_com = CommercialDnaOffer(
    name="Oligo vendor",
    sequence_constraints=[SequenceLengthConstraint(max_length=200)],
    pricing=PerBasepairPricing(0.10),
    lead_time=7,
)

oligo_assembly_station = DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=OligoAssemblyMethod(
        overhang_selector=TmSegmentSelector(
            min_size=15, max_size=25, min_tm=50, max_tm=70
        ),
        min_segment_length=40,
        max_segment_length=200,
        sequence_constraints=[SequenceLengthConstraint(max_length=1500)],
        duration=8,
        cost=2,
    ),
    supplier=oligo_com,
    coarse_grain=20,
    fine_grain=False,
    a_star_factor="auto",
)

gibson_blocks_assembly_station = DnaAssemblyStation(
    name="Gibson Blocks Assembly",
    assembly_method=GibsonAssemblyMethod(
        overhang_selector=FixedSizeSegmentSelector(80),
        min_segment_length=1000,
        max_segment_length=4000,
        duration=8,
        cost=16,
    ),
    supplier=oligo_assembly_station,
    coarse_grain=300,
    fine_grain=False,
    memoize=True,
    a_star_factor="auto",
)

chunks_assembly_station = DnaAssemblyStation(
    name="Chunks assembly (Yeast)",
    assembly_method=GibsonAssemblyMethod(
        overhang_selector=FixedSizeSegmentSelector(300),
        min_segment_length=7000,
        max_segment_length=15000,
        duration=8,
    ),
    supplier=gibson_blocks_assembly_station,
    coarse_grain=1000,
    fine_grain=None,
    logger="bar",
    a_star_factor="auto",
    memoize=True,
)

with open("50kb_sequence.txt", "r") as f:
    sequence = f.read()

print("Generating an assembly plan...")
chunks_assembly_station.prepare_network_on_sequence(sequence)
quote = chunks_assembly_station.get_quote(sequence, with_assembly_plan=True)

print(quote.assembly_step_summary())

print("Generating report...")
assembly_plan_report = quote.to_assembly_plan_report()
assembly_plan_report.write_full_report("report")

print("Done! (see 'report' folder)")

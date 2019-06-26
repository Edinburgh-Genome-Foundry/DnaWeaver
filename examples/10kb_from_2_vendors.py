import dnaweaver as dw

cheap_dna_offer = dw.CommercialDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[
        dw.NoPatternConstraint(enzyme="BsaI"),
        dw.SequenceLengthConstraint(max_length=4000)
    ],
    pricing=dw.PerBasepairPricing(0.10),
)

deluxe_dna_offer = dw.CommercialDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[dw.SequenceLengthConstraint(max_length=3000)],
    pricing=dw.PerBasepairPricing(0.20),
)

assembly_station = dw.DnaAssemblyStation(
    name="Gibson Assembly Station",
    assembly_method=dw.GibsonAssemblyMethod(
        overhang_selector=dw.TmSegmentSelector(min_tm=55, max_tm=70),
        min_segment_length=500,
        max_segment_length=4000
    ),
    supplier=[cheap_dna_offer, deluxe_dna_offer],
    logger='bar',
    coarse_grain=20,
    fine_grain=1
)

sequence = dw.random_dna_sequence(10000, seed=123)
quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
print(quote.assembly_step_summary())
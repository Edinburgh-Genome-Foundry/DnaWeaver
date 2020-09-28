import dnaweaver as dw
import time

cheap_dna_offer = dw.CommercialDnaOffer(
    name="CheapDNA.",
    sequence_constraints=[
        dw.NoPatternConstraint(enzyme="BsaI"),
        dw.SequenceLengthConstraint(max_length=4000),
    ],
    pricing=dw.PerBasepairPricing(0.10),
)

oligo_dna_offer = dw.CommercialDnaOffer(
    name="OliGoo",
    sequence_constraints=[
        dw.GcContentConstraint(min_gc=0.3, max_gc=0.7),
        dw.SequenceLengthConstraint(max_length=100),
    ],
    pricing=dw.PerBasepairPricing(0.07),
    memoize=True,
)

oligo_assembly_station = dw.DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=dw.OligoAssemblyMethod(
        overhang_selector=dw.TmSegmentSelector(
            min_size=15, max_size=25, min_tm=50, max_tm=70
        ),
        min_segment_length=40,
        max_segment_length=200,
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=1500)],
        duration=8,
        cost=30,
    ),
    supplier=oligo_dna_offer,
    coarse_grain=20,
    a_star_factor="auto",
    memoize=True,
)

assembly_station = dw.DnaAssemblyStation(
    name="Gibson Assembly Station",
    assembly_method=dw.GibsonAssemblyMethod(
        overhang_selector=dw.TmSegmentSelector(min_tm=55, max_tm=70),
        min_segment_length=500,
        max_segment_length=4000,
    ),
    supplier=[cheap_dna_offer, oligo_assembly_station],
    logger="bar",
    coarse_grain=100,
    fine_grain=10,
    a_star_factor="auto",
)
print("Looking for the best assembly plan...")
t0 = time.time()
sequence = dw.random_dna_sequence(10000, seed=123)
quote = assembly_station.get_quote(sequence, with_assembly_plan=True)

print(quote.assembly_step_summary())
print("Finished in %.01d seconds" % (time.time() - t0))

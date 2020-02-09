import dnaweaver as dw
import matplotlib.pyplot as plt

# 1. DEFINE THE SUPPLY NETWORK (2 VENDORS FEEDING AN ASSEMBLY STATION)

cheap_dna_offer = dw.CommercialDnaOffer(
    name="CheapDNA",
    sequence_constraints=[
        dw.NoPatternConstraint(enzyme="BsmBI"),
        dw.SequenceLengthConstraint(max_length=4000),
    ],
    pricing=dw.PerBasepairPricing(0.10),
)

deluxe_dna_offer = dw.CommercialDnaOffer(
    name="DeluxeDNA",
    sequence_constraints=[dw.SequenceLengthConstraint(max_length=3000)],
    pricing=dw.PerBasepairPricing(0.20),
)

assembly_stations = [
    dw.DnaAssemblyStation(
        name="Gibson with %s" % supplier.name,
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.TmSegmentSelector(min_tm=55, max_tm=70),
            min_segment_length=500,
            max_segment_length=4000
        ),
        supplier=[supplier],
        logger='bar',
        coarse_grain=20,
        fine_grain=1
    )
    for supplier in [cheap_dna_offer, deluxe_dna_offer]
]
comparator = dw.DnaSuppliersComparator(assembly_stations)

# 2. FIND AN ASSEMBLY PLAN AND CREATE A REPORT.

sequence = dw.load_record("sequence_with_bsmbi_sites.fa")
quote = comparator.get_quote(sequence, with_assembly_plan=True)
print(quote.assembly_step_summary())
assembly_report = quote.to_assembly_plan_report()
_, ax = assembly_report.plot_supply_network()
ax.figure.savefig("supply_network_one_of_two_vendors,pdf", bbox_inches="tight")
figure, ax = plt.subplots(1, figsize=(7, 3))
assembly_report.plot_assembly_blocks(
    ax=ax, parts_offset=0.1, plot_top_assembly=False, legend=True
)
figure.savefig("fragments_from_one_of_2_vendors.pdf", bbox_inches="tight")


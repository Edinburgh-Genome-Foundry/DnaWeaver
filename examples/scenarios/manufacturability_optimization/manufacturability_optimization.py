from dnaweaver import (
    CommercialDnaOffer,
    DnaAssemblyStation,
    GibsonAssemblyMethod,
    DnaSuppliersComparator,
    TmSegmentSelector,
    PerBasepairPricing,
    NoPatternConstraint,
    SequenceLengthConstraint,
)
from dnaweaver.utils import OptimizeManufacturability

import dnachisel

import matplotlib.pyplot as plt

with open("sequence_to_optimize.txt", "r") as f:
    sequence = f.read()

deluxe_dna = CommercialDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[SequenceLengthConstraint(max_length=4000)],
    pricing=PerBasepairPricing(0.20),
    lead_time=10,
)

cheap_dna = CommercialDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[
        NoPatternConstraint(enzyme="BsaI"),
        dnachisel.EnforceGCContent(0.3, 0.7, window=60),
    ],
    pricing=PerBasepairPricing(0.10),
    lead_time=15,
)

# BLOCKS TO CHUNKS ASSEMBLY

gibson_blocks_assembly_station = DnaAssemblyStation(
    name="Gibson Blocks Assembly",
    assembly_method=GibsonAssemblyMethod(
        overhang_selector=TmSegmentSelector(),
        min_segment_length=1000,
        max_segment_length=6000,
        duration=8,
        cost=16,
    ),
    supplier=[deluxe_dna, cheap_dna],
    coarse_grain=30,
    fine_grain=False,
    memoize=True,
    a_star_factor="auto",
)

quote_before = gibson_blocks_assembly_station.get_quote(
    sequence, with_assembly_plan=True
)

print("LOCATING PRICE-DRIVING REGIONS AND OPTIMIZING... PLEASE WAIT")

objective = OptimizeManufacturability(gibson_blocks_assembly_station)

problem = dnachisel.DnaOptimizationProblem(
    sequence=sequence,
    constraints=[dnachisel.EnforceTranslation(location=(0, 9999))],
    objectives=[objective]
)

problem.randomization_threshold = 0  # Forces "random search" mode
problem.max_random_iters = 5
problem.optimize()

print("OPTIMIZATION DONE, GENERATING REPORT")

quote_after = gibson_blocks_assembly_station.get_quote(
    problem.sequence, with_assembly_plan=True
)

fig, axes = plt.subplots(2, figsize=(6, 4))
for title, quote, ax in zip(
    ["Before, optimization", "After optimization"],
    [quote_before, quote_after],
    axes,
):
    report = quote.to_assembly_plan_report()
    ax, _ = report.plot_assembly_blocks(
        parts_offset=0.1,
        plot_top_assembly=False,
        legend=True,
        ax=ax,
        legend_offset=0.15,
    )
    ax.set_title(title, loc="left", fontweight="bold")
fig.tight_layout()
fig.savefig("before_after.pdf")

print("DONE! See result in before_after.pdf")

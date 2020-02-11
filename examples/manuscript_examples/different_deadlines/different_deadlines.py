import os
import dnaweaver as dw
from dnaweaver.biotools import gc_content
import matplotlib.pyplot as plt

oligo_com = dw.CommercialDnaOffer(
    name="Oligo.com",
    sequence_constraints=[dw.SequenceLengthConstraint(max_length=200)],
    pricing=dw.PerBasepairPricing(0.10),
    lead_time=7,
)

deluxe_dna_com = dw.CommercialDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[dw.SequenceLengthConstraint(max_length=10000)],
    pricing=dw.PerBasepairPricing(0.25),
    lead_time=7,
)

cheap_dna_com = dw.CommercialDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[
        dw.SequenceLengthConstraint(max_length=4000),
        dw.NoPatternConstraint(enzyme="AarI"),
        dw.NoPatternConstraint(enzyme="BsaI"),
        lambda seq: (0.4 < gc_content(seq) < 0.6),
    ],
    pricing=dw.PerBasepairPricing(0.10),
    lead_time=15,
)

# OLIGOS TO BLOCKS ASSEMBLY

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
        cost=2,
    ),
    supplier=oligo_com,
    coarse_grain=20,
    fine_grain=False,
    a_star_factor="auto",
    memoize=True,
)

# BLOCKS TO CHUNKS ASSEMBLY

blocks_sources_comparator = dw.DnaSuppliersComparator(
    name="bs_comparator",
    suppliers=[oligo_assembly_station, cheap_dna_com, deluxe_dna_com],
    memoize=True,
)

gibson_blocks_assembly_station = dw.DnaAssemblyStation(
    name="Gibson Blocks Assembly",
    assembly_method=dw.GibsonAssemblyMethod(
        overhang_selector=dw.FixedSizeSegmentSelector(80),
        min_segment_length=1000,
        max_segment_length=4000,
        duration=8,
        cost=16,
    ),
    supplier=blocks_sources_comparator,
    coarse_grain=300,
    fine_grain=False,
    memoize=True,
    a_star_factor="auto",
)

goldengate_blocks_assembly_station = dw.DnaAssemblyStation(
    name="Golden Gate Blocks Assembly",
    assembly_method=dw.GoldenGateAssemblyMethod(
        enzyme="BsmBI",
        wildcard_basepair="A",
        min_segment_length=1000,
        max_segment_length=4000,
        duration=5,
        cost=6,
    ),
    supplier=blocks_sources_comparator,
    coarse_grain=400,
    fine_grain=False,
    memoize=True,
    a_star_factor="auto",
)
ECOLI_DB_PATH = os.path.join("..", "..", "data", "ecoli_blast_db", "ecoli")
ecoli_genome = dw.PcrExtractionStation(
    "E. coli Genome (PCR)",
    primers_supplier=oligo_com,
    homology_selector=dw.TmSegmentSelector(
        min_size=18, max_size=22, min_tm=55, max_tm=65
    ),
    blast_database=ECOLI_DB_PATH,
    max_amplicon_length=10000,
    extra_time=3,
    extra_cost=1,
)

# CHUNKS TO MEGACHUNKS ASSEMBLY

chunks_assembly_station = dw.DnaAssemblyStation(
    name="Chunks assembly (Gibson)",
    assembly_method=dw.GibsonAssemblyMethod(
        overhang_selector=dw.FixedSizeSegmentSelector(300),
        min_segment_length=7000,
        max_segment_length=25000,
        duration=8,
    ),
    supplier=dw.DnaSuppliersComparator(
        [
            ecoli_genome,
            goldengate_blocks_assembly_station,
            gibson_blocks_assembly_station,
            deluxe_dna_com,
        ]
    ),
    coarse_grain=1000,
    fine_grain=None,
    memoize=True,
    a_star_factor="auto",
)

with open("50kb_sequence.txt", "r") as f:
    sequence = f.read()

fig, axes = plt.subplots(1, 4, figsize=(16, 3), sharey=True)
chunks_assembly_station.prepare_network_on_sequence(sequence)
for ax, max_lead_time in zip(axes, [20, 25, 30, 35]):
    quote = chunks_assembly_station.get_quote(
        sequence, max_lead_time=max_lead_time - 1, with_assembly_plan=True
    )

    print("Computing plan for lead time under:", max_lead_time)
    report = quote.to_assembly_plan_report(refine_fragments_locations=False)
    report.plot_assembly_blocks(
        ax=ax, parts_offset=0.1, legend=False, plot_top_assembly=False
    )
    ax.set_title(
        "Best plan under %d days\n\n£%d, %d days"
        % (max_lead_time, quote.price, quote.lead_time)
    )
    ax.set_ylim(top=1)


fig.savefig("different_max_lead_times.pdf")
print ("Done! See different_max_lead_times.pdf")
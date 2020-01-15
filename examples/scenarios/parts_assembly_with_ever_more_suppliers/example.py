import os
from Bio import SeqIO
from plotting_helpers import plot_quote, matplotlib_axes_from_gridspec_array
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
    load_record
)


# DEFINE THE SUPPLY NETWORK

# The EMMA collection of mammalian genetic parts 
emma_collection = GoldenGatePartsLibrary("EMMA", parts_dict={
    record.id: str(record.seq)
    for record in SeqIO.parse("emma_parts.fa", "fasta")
}, memoize=True)

# A medium-price vendor who can provide long parts
company_ingen = CommercialDnaOffer(
    name="InGen",
    pricing=PerBasepairPricing(0.14),
    sequence_constraints=[SequenceLengthConstraint(max_length=2400)]
)

# A cheap vendor who can provide small parts < 1kb
company_tdi = CommercialDnaOffer(
    name="TDI",
    pricing=PerBasepairPricing(0.08),
    sequence_constraints=[SequenceLengthConstraint(max_length=1000)]
)

# An oligos vendor (for oligo assembly and )
company_oligo = CommercialDnaOffer(
    name="Oligo.com",
    pricing=FixedCostPricing(5),
    sequence_constraints=[SequenceLengthConstraint(max_length=100)]
)



mouse_pcr_station = PcrExtractionStation(
    name="E. coli",
    extra_cost=10,
    homology_selector=TmSegmentSelector(),
    primers_supplier=company_oligo,
    blast_database=os.path.join("..", "..", "data", "ecoli_blast_db", "ecoli"),
    memoize=True
)

all_suppliers = [company_ingen,
                 emma_collection,
                 mouse_pcr_station,
                 company_tdi]

assembly_station = DnaAssemblyStation(
    name='GoldenGate Assembly Station',
    assembly_method=GoldenGateAssemblyMethod(
        min_segment_length=40,
        max_segment_length=5000,
        enzyme='BsmBI'
    ),
    supplier=all_suppliers,
    coarse_grain=100,
    logger='bar',
)

for station, color in [(company_ingen, '#fdfdb3'),
                       (company_oligo, '#b8b6e6'),
                       (company_tdi, '#c8f8c4'),
                       (mouse_pcr_station, '#eff5f5'),
                       (emma_collection, '#f8c3c3')]:
    station.report_color = color


# LOAD THE SEQUENCE AND PREPARE THE NETWORK FOR IT. HERE THE
# PREPARATION WILL ONLY MAKE THE PCR STATION PRE-BLAST THE SEQUENCE

record = load_record("example_sequence.gb")
sequence = str(record.seq)
assembly_station.prepare_network_on_sequence(sequence)


# PLOT SUCCESSIVE FIGURES AS THE PROBLEM GETS RESOLVED WITH EVER MORE
# POSSIBLE SUPPLIERS

fig, grid_axes = matplotlib_axes_from_gridspec_array([
    [0, 1, 2, 3],
    [0, 1, 2, 3],
    [4, 5, 6, 7],
], figsize=(14, 6))


for i in range(4):
    # ADD 1 THEN 2... THEN 4 SUPPLIERS
    assembly_station.set_suppliers(all_suppliers[:i+1])
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
    fig, axes = plot_quote(quote, ylim=(-1.2, 1),
                           axes=(grid_axes[i], grid_axes[i + 4]))
fig.tight_layout()
fig.subplots_adjust(hspace=0.2)
fig.savefig("assembly_plans.png", bbox_inches='tight')
print ("Success! See figure assembly_plans.png.")
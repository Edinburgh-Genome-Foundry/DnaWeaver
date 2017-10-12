import pandas as pd
from dnaweaver import (GoldenGateAssemblyMethod, plot_ordering_tree,
                        PartsLibrary, PcrOutStation,
                        DnaAssemblyStation,
                        CommercialDnaOffer, DnaSourcesComparator,
                        random_dna_sequence)
import time

parts = pd.read_csv("all_yeastfab_parts_with_overhangs.csv")
parts_by_sequence = dict(zip(parts.Sequence, parts.Name))
parts_by_name = dict(zip(parts.Name, parts.Sequence))
yeastfab_library = PartsLibrary("YeastFab", parts_dict=parts_by_sequence)


dna_company = CommercialDnaOffer(
    "DNA Company",
    sequence_constraints=[lambda seq: len(seq) < 300],
    price_function=lambda seq: 0.1 * len(seq)
)


yeast_genome = PcrOutStation(
    "Yeast Genome (PCR)",
    primers_dna_source=dna_company,
    blast_database="./cerevisiae_genome/cerevisiae_genome_blast",
    max_amplicon_length=5000
)

gg_station = DnaAssemblyStation(
    "Golden Gate",
    GoldenGateAssemblyMethod(
        "[BsmBI]G",
        min_segment_length=10,
        max_segment_length=3000,
    ),
    dna_source=DnaSourcesComparator([
        yeastfab_library,
        dna_company,
        yeast_genome
    ]),
    a_star_factor=0.02,
    progress_bars=True
)

promoter = parts_by_name['YOR314W_P'][7:-7]
orf = parts_by_name['YOR314W'][11:-7]
terminator = parts_by_name['YOR314W_T'][11:-7]
expression_unit_1 = promoter + orf + terminator
inter_region = random_dna_sequence(800)
with open("./cerevisiae_genome/cerevisiae_genome.fa", "r") as f:
    expression_unit_2 = f.read()[50000:52000].replace("\n", "").upper()
sequence = expression_unit_1 + inter_region + expression_unit_2


t0 = time.time()
yeast_genome.pre_blast(sequence)
quote = gg_station.get_quote(sequence, with_ordering_plan=True)
print quote.ordering_plan.summary()
print time.time() - t0
tree = quote.compute_assembly_tree()

from bokeh.io import output_file
from bokeh.plotting import show
plot = plot_ordering_tree(tree)
output_file("test.html")
show(plot)

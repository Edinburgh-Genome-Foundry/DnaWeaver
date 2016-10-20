from dnaweaver import (ExternalDnaOffer,
                       DnaAssemblyStation,
                       DnaSourcesComparator,
                       PcrOutStation,
                       BuildAGenomeAssemblyMethod,
                       GibsonAssemblyMethod,
                       GoldenGateAssemblyMethod,
                       no_pattern_constraint,
                       random_dna_sequence,
                       gc_content)
from copy import copy
import os


a_star_factor = 0
show_progress = False
ecoli_blast_path = os.path.join("examples_data", "ecoli_blast_db", "ecoli")
sequence_path = os.path.join("examples_data", "multistep_assembly_seq.txt")

oligo_com = ExternalDnaOffer(
    name="Oligo.com",
    sequence_constraints=[lambda seq: len(seq) < 200],
    price_function=lambda sequence: 0.10 * len(sequence),
    lead_time=7
)

deluxe_dna_com = ExternalDnaOffer(
    name="DeluxeDna.com",
    sequence_constraints=[lambda seq: len(seq) < 4000],
    price_function=lambda sequence: 0.20 * len(sequence),
    lead_time=10
)

cheap_dna_com = ExternalDnaOffer(
    name="CheapDna.com",
    sequence_constraints=[
        lambda seq: len(seq) < 4000,
        lambda seq: no_pattern_constraint("GGTCTC"),
        lambda seq: no_pattern_constraint("CACCTGC"),
        lambda seq: (0.4 < gc_content(seq) < 0.6)
    ],
    price_function=lambda sequence: 0.10 * len(sequence),
    lead_time=15,
    memoize=True
)

# OLIGOS TO BLOCKS ASSEMBLY

oligo_assembly_station = DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=BuildAGenomeAssemblyMethod(
        homology_arm_length=20,
        min_segment_length=40,
        max_segment_length=200,
        sequence_constraints=[lambda seq: len(seq) < 1500],
        duration=8,
        cost=5
    ),
    dna_source=oligo_com,
    nucleotide_resolution=20,
    refine_resolution=False,
    progress_bars=show_progress,
    a_star_factor=a_star_factor
)

# BLOCKS TO CHUNKS ASSEMBLY

blocks_sources_comparator = DnaSourcesComparator(
    [
        oligo_assembly_station,
        cheap_dna_com,
        deluxe_dna_com
    ],
    memoize=True
)

gibson_blocks_assembly_station = DnaAssemblyStation(
    name="Gibson Blocks Assembly",
    assembly_method=GibsonAssemblyMethod(
        homology_arm_length=80,
        min_segment_length=1000,
        max_segment_length=4000,
        duration=8,
        cost=8
    ),
    dna_source=blocks_sources_comparator,
    nucleotide_resolution=300,
    refine_resolution=False,
    memoize=True,
    progress_bars=show_progress,
    a_star_factor=a_star_factor
)

goldengate_blocks_assembly_station = DnaAssemblyStation(
    name="Golden Gate Blocks Assembly",
    assembly_method=GoldenGateAssemblyMethod(
        enzyme='BsmBI',
        wildcard_basepair="A",
        min_segment_length=1000,
        max_segment_length=4000,
        avoid_enzyme_in_segments=True,
        duration=5,
        cost=6
    ),
    dna_source=blocks_sources_comparator,
    nucleotide_resolution=400,
    refine_resolution=False,
    memoize=True,
    progress_bars=show_progress,
    a_star_factor=a_star_factor
)

ecoli_genome = PcrOutStation(
    "E. coli Genome (PCR)",
    primers_dna_source=copy(oligo_com),
    blast_database=ecoli_blast_path,
    max_amplicon_length=10000,
    extra_time=3,
    extra_cost=1
)

# CHUNKS TO MEGACHUNKS ASSEMBLY

chunks_assembly_station = DnaAssemblyStation(
    name="Chunks assembly (Gibson)",
    assembly_method=GibsonAssemblyMethod(
        homology_arm_length=300,
        min_segment_length=7000,
        max_segment_length=25000,
        duration=8
    ),
    dna_source=DnaSourcesComparator([
        ecoli_genome,
        gibson_blocks_assembly_station,
        goldengate_blocks_assembly_station
    ]),
    nucleotide_resolution=1000,
    refine_resolution=False,
    progress_bars=show_progress,
    a_star_factor=a_star_factor
)

with open(sequence_path, "r") as f:
    sequence = f.read()
ecoli_genome.pre_blast(sequence)
import time
t0 = time.time()
quote = chunks_assembly_station.get_quote(sequence, with_assembly_plan=True)
print(quote.assembly_step_summary())
print("Finished in %.02fs" % (time.time()-t0))

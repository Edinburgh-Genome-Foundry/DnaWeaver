from dnaweaver import (DnaOffer, DnaAssemblyOffer, DnaOrderingProblem,
                        GoldenGateAssemblyMethod, GibsonAssemblyMethod)
from dnachisel import random_dna_sequence
import dnachisel.constraints as cst
import numpy


# OLIGO COMPANIES

big_dna_com = DnaOffer(
    name="BigDNA.com",
    pricing=lambda sequence: 0.25 * len(sequence),
    constraints=[lambda seq: len(seq) < 4000], memoize=True
)


# BLOCKS ASSEMBLY

blocks_assembler_gibson_assembly_offer = DnaAssemblyOffer(
    "Blocks Assembler",
    DnaOrderingProblem(offers=[big_dna_com],
                       assembly_method=GibsonAssemblyMethod(40)),
    memoize=True,
    min_segment_length=2000,
    max_segment_length=4000,
    nucleotide_resolution=200,
    refine_resolution=False
)

blocks_assembler_golden_gate_offer = DnaAssemblyOffer(
    "Blocks Assembler",
    DnaOrderingProblem(offers=[big_dna_com],
                       assembly_method=GoldenGateAssemblyMethod()),
    memoize=True,
    min_segment_length=2000,
    max_segment_length=4000,
    nucleotide_resolution=200,
    refine_resolution=False
)


# CHUNKS ASSEMBLY

chunks_assembler_offer = DnaAssemblyOffer(
    "Chunks Assembler",
    DnaOrderingProblem(offers=[blocks_assembler_gibson_assembly_offer,
                               blocks_assembler_golden_gate_offer],
                       assembly_method=GibsonAssemblyMethod(300)),
    memoize=True,
    min_segment_length=15000,
    max_segment_length=25000,
    nucleotide_resolution=1000,
    refine_resolution=False
)

numpy.random.seed(1234)
sequence = random_dna_sequence(50000)
chunks_assembler_offer.best_ordering_plan(sequence)

# WRITE GENBANK FILE

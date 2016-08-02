from dnaweaver import (DnaOffer, DnaAssemblyOffer, DnaOrderingProblem,
                        BuildAGenomeAssemblyMethod, GibsonAssemblyMethod)
from dnachisel import random_dna_sequence
import dnachisel.constraints as cst
import numpy


# OLIGO COMPANIES

cheap_dna_com = DnaOffer(
    name="CheapDNA.com",
    constraints=[cst.GCContentConstraint(0.4, 0.6, gc_window=50)],
    pricing=lambda sequence: 0.10 * len(sequence),
    memoize=True
)

deluxe_dna_com = DnaOffer(
    name="DeluxeDNA.com",
    pricing=lambda sequence: 0.20 * len(sequence),
    constraints=[], memoize=True
)


# OLIGO ASSEMBLY

oligos_assembler_offer = DnaAssemblyOffer(
    "Oligos Assembler",
    DnaOrderingProblem(offers=[cheap_dna_com, deluxe_dna_com],
                       assembly_method=BuildAGenomeAssemblyMethod(20)),
    memoize=True,
    min_segment_length=40,
    max_segment_length=100,
    nucleotide_resolution=10,
    refine_resolution=False
)


# BLOCKS ASSEMBLY

blocks_assembler_offer = DnaAssemblyOffer(
    "Blocks Assembler",
    DnaOrderingProblem(offers=[oligos_assembler_offer],
                       assembly_method=GibsonAssemblyMethod(40)),
    memoize=True,
    min_segment_length=2000,
    max_segment_length=4000,
    nucleotide_resolution=200,
    refine_resolution=False
)


# CHUNKS ASSEMBLY

chunks_assembler_offer = DnaAssemblyOffer(
    "Chunks Assembler",
    DnaOrderingProblem(offers=[blocks_assembler_offer],
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

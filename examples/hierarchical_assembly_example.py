"""Hierarchical assembly.

Problem
--------


What we show
------------

"""

from dnaadvisor import (DnaOffer, DnaOrderingProblem,
                        BuildAGenomeAssemblyMethod, GibsonAssemblyMethod)
from dnachisel import random_dna_sequence
import dnachisel.constraints as cst
import primer3   # primer3 Python binding
import numpy

numpy.random.seed(1234)
sequence = random_dna_sequence(50000)


# OLIGO ASSEMBLY

def has_melting_temperature_between_50_and_55(location):
    """Return False if the 20-basepair segment around the location
    has a melting temperature outside 50-55 Celsius."""
    if min(location, len(sequence) - location) < 20:
        return True
    subsequence = sequence[location - 10:location + 10]
    melting_temperature = primer3.calcTm(subsequence)
    return 40 < melting_temperature < 60

cheap_dna_com = DnaOffer(name="CheapDNA.com",
                         constraints=[cst.GCContentConstraint(
                             0.4, 0.6, gc_window=50)],
                         pricing=lambda sequence: 0.10 * len(sequence),
                         memoize=True)

deluxe_dna_com = DnaOffer(name="DeluxeDNA.com", constraints=[],
                          pricing=lambda sequence: 0.20 * len(sequence),
                          memoize=True)

oligos_assembler_problem = DnaOrderingProblem(
    offers=[cheap_dna_com, deluxe_dna_com],
    location_filters=(has_melting_temperature_between_50_and_55,),
    assembly_method=BuildAGenomeAssemblyMethod(20)
)
oligos_assembler_offer = oligos_assembler_problem.as_DnaOffer(
    "Oligos Assembler",
    memoize=True,
    min_segment_length=40,
    max_segment_length=100,
    nucleotide_resolution=10,
    refine_resolution=False
)


# BLOCKS ASSEMBLY

blocks_assembler_problem = DnaOrderingProblem(
    offers=[oligos_assembler_offer],
    assembly_method=GibsonAssemblyMethod(40),
)
blocks_assembler_offer = blocks_assembler_problem.as_DnaOffer(
    "Blocks Assembler",
    memoize=True,
    min_segment_length=2000,
    max_segment_length=4000,
    nucleotide_resolution=200,
    refine_resolution=False
)


# CHUNKS ASSEMBLY

chunks_assembly_problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[blocks_assembler_offer],
    assembly_method=GibsonAssemblyMethod(300)
)


solution = chunks_assembly_problem.solve(
    min_segment_length=15000, max_segment_length=25000,
    nucleotide_resolution=1000,
    refine_resolution=False
)

# WRITE GENBANK FILE

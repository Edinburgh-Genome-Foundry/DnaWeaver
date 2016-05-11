"""Example of oligo design for oligo assembly using DnaAdvisor

In this example we will order DNA oligos that will be assembled using an
oligo assembly method
"""

from dnaadvisor import DnaOffer, DnaOrderingProblem, BuildAGenomeAssemblyMethod
from dnachisel import random_dna_sequence
import primer3   # primer3 Python binding
import RNA       # ViennaRNA binding
import numpy

numpy.random.seed(1234)

sequence = random_dna_sequence(1500)
assembly_method = BuildAGenomeAssemblyMethod(20)


def has_melting_temperature_between_50_and_55(location):
    """Return False if the 20-basepair segment around the location
    has a melting temperature outside 50-55 Celsius."""
    if min(location, len(sequence) - location) < 20:
        return True
    subsequence = sequence[location - 10:location + 10]
    melting_temperature = primer3.calcTm(subsequence)
    return 50 < melting_temperature < 55


def has_weak_secondary_structure(segment):
    """Return False if the segment corresponds to an oligo with a
    secondary structure that is too stable"""
    fragment = assembly_method.compute_fragment_sequence(segment, sequence)
    return RNA.fold(fragment)[1] > -40


company_1 = DnaOffer(name="Company 1", constraints=[],
                     pricing=lambda sequence: 0.10 * len(sequence))

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[company_1],
    location_filters=(has_melting_temperature_between_50_and_55,),
    segment_filters=(has_weak_secondary_structure,),
    assembly_method=assembly_method
)


solution = problem.solve(
    min_segment_length=40,
    max_segment_length=100,
    nucleotide_resolution=4,
)

print (solution.summary())

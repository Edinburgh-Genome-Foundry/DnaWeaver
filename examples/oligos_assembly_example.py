"""Example of oligo design for oligo assembly using DnaAdvisor

Problem
--------
We consider a random 1500bp DNA sequence that we wish to assemble using an
oligo assembly method, e.g. Build-a-Genome. The oligos should verify:
- Length between 40 and 100
- 20bp overlaps between adjacent oligos
- Overlaps should have a Tm between 50-55 degC (as measured by primer3)
- Oligos should have a weak secondary structure (dG > -40 as measured
  by ViennaRNA)

What we show
------------
- We show how to use `locations_filter` to define the cut locations that will
  produce overlaps with the right Tm
- We show how to use `segments_filter` to filter out segments which will
  produce oligos with too strong a secondary structure.
- We show that DnaAdvisor can come to a nice solution under these conditions.

Technical note
--------------
Note that edges computing (in particular the secondary structure with
ViennaRNA) is the bottleneck of computations.
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


oligo_company = DnaOffer(name="Oligos Company", constraints=[],
                         pricing=lambda sequence: 0.10 * len(sequence))

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[oligo_company],
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

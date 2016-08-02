"""Example of oligo design for oligo assembly using DnaWeaver

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
- We show that DnaWeaver can come to a nice solution under these conditions.

Technical note
--------------
Note that edges computing (in particular the secondary structure with
ViennaRNA) is the bottleneck of computations.
"""


from dnaweaver import (ExternalDnaOffer, DnaAssemblyStation,
                        BuildAGenomeAssemblyMethod)
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


assembly_station = DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=BuildAGenomeAssemblyMethod(
        homology_arm_length=20,
        min_segment_length=40,
        max_segment_length=100,
        duration=7
    ),
    dna_source=ExternalDnaOffer(
        name="Oligos Company",
        pricing=lambda seq: 0.10 * len(seq)
    ),
    nucleotide_resolution=10,
    refine_resolution=1
)

quote = assembly_station.get_quote(sequence, with_ordering_plan=True)

print (quote)
if quote.accepted:
    print (quote.ordering_plan.summary())

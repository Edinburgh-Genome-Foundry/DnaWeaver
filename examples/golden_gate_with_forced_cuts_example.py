"""Advanced Golden Gate cutting example for DnaWeaver.

This example is the continuation of golden_gate_simple_examples.py.

In this example we consider a random 5000bp DNA sequence that we wish to order
to the company InGen.

The problem is that our sequence has an AarI restriction site as well as
an homopolymer of nine "G"s, which InGen forbi
As a solution, we locate these sites and forced the solver to cut in the middle
of these sites.
"""

from dnaweaver import (CommercialDnaOffer, DnaAssemblyStation,
                        GoldenGateAssemblyMethod)
from dnachisel import (random_dna_sequence, enzyme_pattern,
                       NoPatternConstraint, homopolymer_pattern)
import numpy

numpy.random.seed(1234)

sequence = random_dna_sequence(5000)

forbidden_patterns = [
    enzyme_pattern("AarI"),
    homopolymer_pattern("G", 9)
]

for pattern in forbidden_patterns:
    matches = pattern.find_matches(sequence)
    print ("Pattern %s found at locations %s" % (pattern, matches))

def force_cuts(sequence):
    all_forbidden_patterns_centers = sorted([
        (a + b) // 2
        for pattern in forbidden_patterns
        for (a, b) in pattern.find_matches(sequence)
    ])
    return all_forbidden_patterns_centers

all_forbidden_patterns_centers = sorted([
    (a + b) // 2
    for pattern in forbidden_patterns
    for (a, b) in pattern.find_matches(sequence)
])

cheap_dna_com = CommercialDnaOffer(
    "CheapDNA.com",
    sequence_constraints= [NoPatternConstraint(pattern)
                           for pattern in forbidden_patterns],
    price_function=lambda sequence: 0.10 * len(sequence))


assembly_station = DnaAssemblyStation(
    "Golden Gate Assembly Station",
    assembly_method=GoldenGateAssemblyMethod(
        left_overhang="[BsaI]A",
        min_segment_length=100,
        max_segment_length=3500,
        force_cuts=force_cuts
    ),
    dna_source=cheap_dna_com,
    nucleotide_resolution=10,
    refine_resolution=1
)

quote = assembly_station.get_quote(sequence, with_ordering_plan=True)

print (quote)
if quote.accepted:
    print (quote.ordering_plan.summary())

# This will print:
# ----------------
#
# Pattern [C][A][C][C][T][G][C] (AarI) found at locations [(1237, 1244)]
# Pattern [G][G][G][G][G][G][G][G][G] found at locations [(3247, 3256)]
# Ordering plan:
#   (0, 1240) Company 1 127.00$
#   (1240, 3251) Company 1 204.30$
#   (3251, 5000) Company 1 177.90$
#   Total:509$

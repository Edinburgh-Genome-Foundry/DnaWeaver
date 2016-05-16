"""Advanced Golden Gate cutting example for DnaAdvisor.

This example is the continuation of golden_gate_simple_examples.py.

In this example we consider a random 5000bp DNA sequence that we wish to order
to the company InGen.

The problem is that our sequence has an AarI restriction site as well as
an homopolymer of nine "G"s, which InGen forbids.
As a solution, we locate these sites and forced the solver to cut in the middle
of these sites.
"""

from dnaadvisor import *
from dnachisel import (random_dna_sequence, enzyme_pattern,
                       NoPatternConstraint, homopolymer_pattern)
from dnachisel.biotools import gc_content
import numpy as np

np.random.seed(1234)
sequence = random_dna_sequence(5000)

forbidden_patterns = [
    enzyme_pattern("AarI"),
    homopolymer_pattern("G", 9)
]

for pattern in forbidden_patterns:
    matches = pattern.find_matches(sequence)
    print ("Pattern %s found at locations %s" % (pattern, matches))

all_forbidden_patterns_centers = sorted([
    (a + b) // 2
    for pattern in forbidden_patterns
    for (a, b) in pattern.find_matches(sequence)
])

company_1 = DnaOffer(name="Company 1",
                     constraints=[NoPatternConstraint(pattern)
                                  for pattern in forbidden_patterns],
                     pricing=lambda sequence: 0.10 * len(sequence))

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[company_1],
    assembly_method=GoldenGateAssemblyMethod()
)

solution = problem.solve(
    min_segment_length=100,
    max_segment_length=3500,
    nucleotide_resolution=50,
    forced_cuts=all_forbidden_patterns_centers,
    refine_resolution=1
)

print solution.summary()

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

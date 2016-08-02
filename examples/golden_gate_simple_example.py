"""Simple Golden Gate cutting example for DnaWeaver.

Problem
--------
we consider a random 5000bp DNA sequence that we wish to order to the company
InGen.

The problem is that our sequence has an AarI restriction site, which InGen
forbids, around position 1240. The solution is to order 2 segments of DNA,
[0,1240] and [1240,5000], so that none of the segments feature the site, then
the segments are assembled by Gibson Assembly.

What we show
------------
- Here we show that DnaWeaver can come to that solution with no guidance.
- We also demonstrate the use of `locations_filters` to filter out locations
  which would be unsuitable for Golden Gate assembly cuts (with a GC content of
  100%).
- The last line also shows how to print the solution to a spreadsheet, ready to
  order !

Technical note
--------------
If we increase the nucleotide resolution of our solver to 50 (to
make it faster) the solver will overlook the cut in location 1240, and
therefore come to the conclusion that there is no solution to the problem.
In and other example (golden_gate_with_forced_cuts_example.py) we show a
smarter way to solve this problem by forcing the location of certain cuts.
"""

from dnaweaver import *
from dnachisel import random_dna_sequence, enzyme_pattern, NoPatternConstraint
from dnachisel.biotools import gc_content
import numpy as np

np.random.seed(1234)

sequence = random_dna_sequence(5000)

aarI_site = enzyme_pattern("AarI")
sites_locations = aarI_site.find_matches(sequence)
enzyme_sites_centres = sorted([(a + b) // 2 for (a, b) in sites_locations])
print ("BsaI site found around positions %s" % sites_locations)

company_1 = DnaOffer(
    name="Company InGen",
    constraints=[NoPatternConstraint(aarI_site)],
    pricing=lambda sequence: 0.10 * len(sequence)
)

def is_suitable_gg_linker(location):
    """Return True iff the 4bp site around `location` has less than 100% GC"""
    subsequence = sequence[location - 2:location + 2]
    return (len(subsequence) == 4) and (gc_content(subsequence) < 1)

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[company_1],
    location_filters=(is_suitable_gg_linker,),
    assembly_method=GoldenGateAssemblyMethod()
)

solution = problem.solve(
    min_segment_length=100,
    max_segment_length=4000,
    nucleotide_resolution=10,
)

print (solution.summary())

# Save as csv
solution.to_dataframe().to_csv("golden_gate_ordering_sheet.csv", index=False)

# This will print:
# ----------------
#
# BsaI site found around positions [(1237, 1244)]
# Ordering plan:
#   (0, 1240) Company InGen 125.60$
#   (1240, 5000) Company InGen 377.60$
#   Total:503$

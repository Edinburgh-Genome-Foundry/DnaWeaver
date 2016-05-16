"""Basic example for Dna Advisor: ordering DNA from 2 different companies

Problem
--------
We consider a random 5000bp DNA sequence that we wish to assmeble using
Gibson Assembly. The fragments, which will have 40bp overlap, can be ordered
from to two different companies:
- Company 1 charges 0.1$/bp for any sequence up to 4000bp without BsaI site
- Company 2 charges 0.2$/bp for any sequence up to 3000bp

As our sequence has 2 BsaI sites, some fragments will have to be ordered to
Company 2, but they should be as short as possible as they are more expensive.
The obvious solution is to use company 2 for just the fragments around BsaI
sites, and company 1 for the rest.

What we show
------------
We show that DnaAdvisor can come to the best solution with no guidance.

Technical note
--------------
Try playing with the `cuts_number_penalty` parameter in the problem definition
to reduce the number of cuts in the proposed solution.
For instance if you increase `cuts_number_penalty` to 500 you will see the
number of segments to order fall from 7 to just 3.
"""

from dnaadvisor import *
from dnachisel import (random_dna_sequence, enzyme_pattern,
                       NoPatternConstraint, SequenceLengthConstraint)
import numpy as np

np.random.seed(123)  # Set the randomizer's seed to ensure reproducibility.

sequence = random_dna_sequence(10000)

enzyme_site = enzyme_pattern("BsaI")
print ("BsaI site found at positions %s" % enzyme_site.find_matches(sequence))

company_1 = DnaOffer(
    name="Company 1",
    constraints=[NoPatternConstraint(enzyme_site),
                 SequenceLengthConstraint(max_length=4000)],
    pricing=lambda sequence: 0.10 * len(sequence)
)

company_2 = DnaOffer(
    name="Company 2",
    constraints=[SequenceLengthConstraint(max_length=3000)],
    pricing=lambda sequence: 0.20 * len(sequence)
)

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[company_1, company_2],
    assembly_method=GibsonAssemblyMethod(40),
    cuts_number_penalty=500
)

solution = problem.solve(
    min_segment_length=100,
    max_segment_length=4000,
    nucleotide_resolution=20,
    refine_resolution=1,

)

print (solution.summary())

# This will print:
# ----------------
#
# BsaI site found at positions [(4435, 4441), (5307, 5313)]
# Ordering plan:
#   (0, 449) Company 1 48.90$
#   (449, 4369) Company 1 400.00$
#   (4369, 4476) Company 2 37.40$
#   (4476, 5249) Company 1 85.30$
#   (5249, 5350) Company 2 36.20$
#   (5350, 6040) Company 1 77.00$
#   (6040, 10000) Company 1 400.00$
#   Total:1084$

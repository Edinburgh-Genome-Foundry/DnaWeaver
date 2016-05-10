"""
Apply cuts_number_penalty=300
"""

from dnaadvisor import *
from dnachisel import random_dna_sequence, enzyme_pattern, NoPatternConstraint
import numpy as np

np.random.seed(123)

sequence = random_dna_sequence(10000)

enzyme_site = enzyme_pattern("BsaI")


company_1 = DnaOffer(
    name="Company 1",
    constraints=[NoPatternConstraint(enzyme_site)],
    pricing=lambda sequence: 0.10 * len(sequence)
)

company_2 = DnaOffer(
    name="Company 2",
    constraints=[],
    pricing=lambda sequence: 0.20 * len(sequence)
)

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[company_1, company_2],
    assembly_method=GibsonAssemblyMethod(20)
)

best_cut = optimize_costs_with_graph(problem, segment_length_range=(150, 4000),
                                     nucleotide_resolution=100, refine=False,
                                     cuts_number_penalty=0)

print ("BsaI site found at positions %s" % enzyme_site.find_matches(sequence))
print ("Here is the proposed ordering solution:")
offers = sorted(best_cut.values(), key=lambda o: o.zone)
for offer in  offers:
    print (offer)
print "Total price: %d $" % sum(o.price for o in offers)

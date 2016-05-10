"""
Apply cuts_number_penalty=300
"""

from dnaadvisor import *
from dnachisel import random_dna_sequence, SequenceLengthConstraint
import numpy as np
np.random.seed(123)


def filter_has_moderate_melting_temperature(location, sequence):
    if min(location, len(sequence) - location) < 20:
        return True
    subsequence = sequence[location - 10:location + 10]
    melting_temperature = sum({"A": 2, "T": 2, "G": 4, "C": 4}[base]
                              for base in subsequence)
    return 60 < melting_temperature < 66

company_1 = DnaOffer(
    name="Company 1",
    constraints=[SequenceLengthConstraint(max_length=4000)],
    pricing=lambda sequence: 0.10 * len(sequence)
)

problem = DnaOrderingProblem(
    sequence=random_dna_sequence(10000),
    offers=[company_1],
    assembly_method=GibsonAssemblyMethod(20),
)

best_cut = optimize_costs_with_graph(
    problem, segment_length_range=(150, 4000),
    location_filter=filter_has_moderate_melting_temperature,
    nucleotide_resolution=100, refine=False
)

print ("Proposed ordering solution:")
offers = sorted(best_cut.values(), key=lambda o: o.zone)
for offer in offers:
    print (offer)
print "Total price: %d $" % sum(o.price for o in offers)

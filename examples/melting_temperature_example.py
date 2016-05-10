"""
Apply cuts_number_penalty=300
"""

from dnaadvisor import *
from dnachisel import random_dna_sequence, SequenceLengthConstraint
import numpy as np
np.random.seed(123)

sequence = random_dna_sequence(10000)

def valid_melting_temperature(location):
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
    sequence= sequence,
    offers=[company_1],
    location_filters = (valid_melting_temperature,),
    assembly_method=GibsonAssemblyMethod(20),
)


solution = problem.solve(
    min_segment_length=150,
    max_segment_length=4000,
    nucleotide_resolution=100,
    refine_resolution=False
)

print (problem.ordering_plan_summary(solution))

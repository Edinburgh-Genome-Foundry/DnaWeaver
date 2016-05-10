from dnaadvisor import *
from dnachisel import random_dna_sequence, PROVIDERS_CONSTRAINTS
import numpy as np

np.random.seed(123)

sequence = (random_dna_sequence(200) + 'GGTCTC' +
            random_dna_sequence(300) + 'GGTCTC' +
            random_dna_sequence(200))

gen9_offer = DnaOffer(
    name="Gen9Offer",
    constraints = PROVIDERS_CONSTRAINTS["Gen9"],
    pricing = lambda sequence: 0.08*len(sequence)
)

idt_offer = DnaOffer(
    name="IDTOffer",
    constraints= PROVIDERS_CONSTRAINTS["IDT"],
    pricing = lambda sequence: 0.10*len(sequence)
)
problem = DnaOrderingProblem(
    sequence= sequence,
    offers = [gen9_offer, idt_offer],
    assembly_method= GibsonAssemblyMethod(20)
)

offers =  optimize_costs_with_graph(problem, segment_length_range=(5, 400),
                                    nucleotide_resolution=10)
offers = sorted(offers.values(), key=lambda o: o.zone)
for offer in offers:
    print (offer)
print "total: %d $" % sum(o.price for o in offers)

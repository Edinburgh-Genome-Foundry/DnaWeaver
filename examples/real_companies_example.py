from dnaadvisor import *
from dnachisel import random_compatible_dna_sequence, PROVIDERS_CONSTRAINTS
import numpy as np

np.random.seed(12345)

gen9_rules = PROVIDERS_CONSTRAINTS["Gen9"]
idt_rules = PROVIDERS_CONSTRAINTS["IDT"]

print ("==> Gen9 rules:\n%s\n" % "\n".join(map(str, gen9_rules)))
print ("==> IDT rules:\n%s\n" % "\n".join(map(str, idt_rules)))

# NOTE: The seed is set so as to make sure that the "random" fragments and
# their subfragments are indeed gen9-compatible
sequence = (random_compatible_dna_sequence(100, gen9_rules) + 'GGTCTC' +
            random_compatible_dna_sequence(500, gen9_rules) + 'GGTCTC' +
            random_compatible_dna_sequence(100, gen9_rules))

gen9_offer = DnaOffer(name="Gen9Offer", constraints=gen9_rules,
                      pricing=lambda sequence: 0.10 * len(sequence))

idt_offer = DnaOffer(name="IDTOffer", constraints=idt_rules,
                     pricing=lambda sequence: 0.20 * len(sequence))

problem = DnaOrderingProblem(
    sequence=sequence,
    offers=[gen9_offer, idt_offer],
    assembly_method=GibsonAssemblyMethod(20)
)

solution = problem.solve(
    min_segment_length=50,
    max_segment_length=700,
    nucleotide_resolution=5,
)

print (solution.summary())

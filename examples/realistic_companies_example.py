"""An example of cutting with DNA Advisor using real DNA companies rules .

Problem
--------

We want to order fragments to assemble a piece of DNA, we consider two
fictitious companies which have many complex rules like in real life:

Company Cell9 charges 0.10$/bp and applies the following constraints:

- GC content must be between 40-65% (and 25-80%, over 50bp windows)
- No BsaI or AarI sites.
- No homopolymers of size 8+ (for A,T,C) or 6+ (for G).

Company TDI charges 0.20$/bp and applies the following constraints:

- GC content must be:
    - 40-68% globally
    - 28-76% over 100bp windows
    - 15-90% over 20bp windows
    - 24-76% in the 30bp terminal windows
- No 3-mers repeated 5+ times, or 2-mers repeated 9+ times
- No 6+ homopolymers of G-C or 9+ homopolymers of A-T
- No hairpins (defined as 20-mers with a reverse-complement in the next 200bp window)
- No homopolymers of size 8+ (for A,T,C) or 6+ (for G).

What we show
------------
- We show that the company constraints can be modelled using DnaChisel
- We show that in this example DnaAdvisor finds the smartest solution.
"""

from dnaadvisor import *
from dnachisel import random_compatible_dna_sequence, patterns
import dnachisel.constraints as cst
import numpy as np

np.random.seed(12345)

cell9_rules = [
    cst.NoPatternConstraint(patterns.enzyme_pattern("BsaI")),
    cst.NoPatternConstraint(patterns.enzyme_pattern("AarI")),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("A", 9)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("T", 9)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("G", 6)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("C", 9)),
    cst.GCContentConstraint(0.4, 0.65),
    cst.GCContentConstraint(0.25, 0.80, gc_window=50)
]
tdi_rules = [
    cst.GCContentConstraint(0.25, 0.68),
    cst.GCContentConstraint(0.28, 0.76, gc_window=100),
    cst.GCContentConstraint(0.15, 0.90, gc_window=20),
    cst.TerminalGCContentConstraint(0.24, 0.76, window_size=30),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("A", 13)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("T", 13)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("G", 6)),
    cst.NoPatternConstraint(patterns.homopolymer_pattern("C", 6)),
    cst.NoPatternConstraint(patterns.repeated_kmers(3, n_repeats=5)),
    cst.NoPatternConstraint(patterns.repeated_kmers(2, n_repeats=9)),
    cst.NoHairpinsIDTConstraint(stem_size=20, hairpin_window=200)
]

cell9_offer = ExternalDnaOffer(
    name="Cell9 offer",
    sequence_constraints=cell9_rules,
    price_function=lambda sequence: 0.10 * len(sequence)
)

tdi_offer = ExternalDnaOffer(
    name="TDI Offer",
    sequence_constraints=tdi_rules,
    price_function=lambda sequence: 0.10 * len(sequence)
)


# NOTE: The seed is set so as to make sure that the "random" fragments and
# their subfragments are indeed gen9-compatible
sequence = (random_compatible_dna_sequence(1000, cell9_rules) + 'GGTCTC' +
            random_compatible_dna_sequence(5000, cell9_rules) + 'GGTCTC' +
            random_compatible_dna_sequence(1000, cell9_rules))

assembly_station = DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=BuildAGenomeAssemblyMethod(
        homology_arm_length=20,
        min_segment_length=50,
        max_segment_length=700,
        duration=7
    ),
    dna_source=DnaSourcesComparator([cell9_offer, tdi_offer]),
    nucleotide_resolution=20,
    refine_resolution=False
)

quote = assembly_station.get_quote(sequence, with_ordering_plan=True)

print (quote)
if quote.accepted:
    print (quote.ordering_plan.summary())

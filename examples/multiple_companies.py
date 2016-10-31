"""Basic example for Dna Weaver: ordering DNA from 2 different companies

Problem
--------
We consider a random 5000bp DNA sequence that we wish to assemble using
Gibson Assembly. The fragments, which will have 40bp overlap, can be ordered
from to two different companies:
- CheapDNA.com charges 0.1$/bp for any sequence up to 4000bp without BsaI site
- DeluxeDNA.com charges 0.2$/bp for any sequence up to 3000bp

As our sequence has 2 BsaI sites, some fragments will have to be ordered to
DeluxeDNA.com, but they should be as short as possible as they are more
expensive.
The obvious solution is to use company 2 for just the fragments around BsaI
sites, and company 1 for the rest.

What we show
------------
We show that DnaWeaver can come to the best solution with no guidance.

Technical note
--------------
Try playing with the `cuts_number_penalty` parameter in the problem definition
to reduce the number of cuts in the proposed solution.
For instance if you increase `cuts_number_penalty` to 500 you will see the
number of segments to order fall from 7 to just 3.
"""
from dnaweaver import *

cheap_dna_offer = ExternalDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[
        NoPatternConstraint("GGTCTC"),
        lambda seq: len(seq) < 4000
    ],
    price_function=PerBasepairPricing(per_basepair_price=0.10),
)

deluxe_dna_offer = ExternalDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[lambda seq: len(seq) < 3000],
    price_function=PerBasepairPricing(per_basepair_price=0.20),
)

assembly_station = DnaAssemblyStation(
    name="Gibson Assembly Station",
    assembly_method=GibsonAssemblyMethod(homology_arm_length=20,
                                         min_segment_length=500,
                                         max_segment_length=4000),
    dna_source=DnaSourcesComparator([cheap_dna_offer, deluxe_dna_offer]),
    nucleotide_resolution=10,
    refine_resolution=False
)


print("Now finding a price-optimal assembly strategy for a 10kb sequence")
sequence = random_dna_sequence(10000, seed=123)
import time
t0 = time.time()
quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
print (quote.assembly_step_summary())
print("Finished in %.02fs" % (time.time()-t0))

# This will print:
# ----------------
#
# Ordering plan:
#   (0, 500): From CheapDNA.com, price 52.00
#   (500, 4420): From CheapDNA.com, price 396.00
#   (4420, 5330): From DeluxeDNA.com, price 190.00
#   (5330, 6030): From CheapDNA.com, price 74.00
#   (6030, 10000): From CheapDNA.com, price 399.00
#   Price:1111.00

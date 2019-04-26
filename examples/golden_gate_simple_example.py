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

from dnaweaver import (CommercialDnaOffer, GoldenGateAssemblyMethod,
                       PerBasepairPricing, DnaAssemblyStation)
from dnachisel import random_dna_sequence, enzyme_pattern, AvoidPattern

company = CommercialDnaOffer(
    name="Company InGen",
    sequence_constraints=[AvoidPattern(enzyme='AarI')],
    pricing=PerBasepairPricing(0.08)
)
assembly_station = DnaAssemblyStation(
    name='GoldenGate Assembly Station',
    assembly_method=GoldenGateAssemblyMethod(
        min_segment_length=50,
        max_segment_length=2000
    ),
    dna_source=company,
    coarse_grain=50,
    logger='bar'
)

sequence = random_dna_sequence(4000, seed=123)
sites_locations = enzyme_pattern("AarI").find_matches(sequence)
enzyme_sites_centres = [(l.start + l.end) // 2 for l in sites_locations]
print ("AarI site(s) were found around position(s) %s" % enzyme_sites_centres)
quote = assembly_station.get_quote(sequence,  with_assembly_plan=True)

print (quote)
if quote.accepted:
    print (quote.assembly_step_summary())

# This will print:
# ----------------
#
# BsaI site found around positions [(1237, 1244)]
# Ordering plan:
#   (0, 1240) Company InGen 125.60$
#   (1240, 5000) Company InGen 377.60$
#   Total:503$

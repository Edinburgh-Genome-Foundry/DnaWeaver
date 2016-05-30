from dnaadvisor import (ExternalDnaOffer,
                        DnaAssemblyStation,
                        GibsonAssemblyMethod,
                        BuildAGenomeAssemblyMethod,
                        DnaSourcesComparator)
from dnachisel import random_dna_sequence
import dnachisel.constraints as cst
import numpy


# OLIGO COMPANIES

cheap_dna_com = ExternalDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[cst.GCContentConstraint(0.4, 0.6, gc_window=50),
                          cst.SequenceLengthConstraint(max_length=200)],
    price_function=lambda sequence: 0.10 * len(sequence),
    lead_time=10,
    memoize=True
)

deluxe_dna_com = ExternalDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[cst.SequenceLengthConstraint(max_length=200)],
    price_function=lambda sequence: 0.20 * len(sequence),
    lead_time=5,
    memoize=True
)

big_dna_com = ExternalDnaOffer(
    name="BigDNA.com",
    sequence_constraints=[cst.SequenceLengthConstraint(2000, 4000)],
    price_function=lambda sequence: 0.40 * len(sequence),
    lead_time=10,
    memoize=True
)


# OLIGOS TO BLOCKS ASSEMBLY

oligo_assembly_station = DnaAssemblyStation(
    name="Oligo Assembly Station",
    assembly_method=BuildAGenomeAssemblyMethod(
        homology_arm_length=20,
        min_segment_length=40,
        max_segment_length=100,
        duration=8,
        cost=0
    ),
    dna_source=DnaSourcesComparator([
        cheap_dna_com,
        deluxe_dna_com
    ]),
    nucleotide_resolution=10,
    refine_resolution=False,
    memoize=True,
    progress_bars=False,
    a_star_factor=0
)


# BLOCKS TO CHUNKS ASSEMBLY

blocks_assembly_station = DnaAssemblyStation(
    name="Blocks Assembly Station",
    assembly_method=GibsonAssemblyMethod(
        homology_arm_length=40,
        min_segment_length=2000,
        max_segment_length=4000,
        duration=8
    ),
    dna_source=DnaSourcesComparator([
        oligo_assembly_station,
        big_dna_com
    ]),
    nucleotide_resolution=200,
    refine_resolution=False,
    memoize=True,
    progress_bars=False,
    a_star_factor=0
)

# CHUNKS TO MEGACHUNKS ASSEMBLY

chunks_assembly_station = DnaAssemblyStation(
    name="Chunks Assembly Station",
    assembly_method=GibsonAssemblyMethod(
        homology_arm_length=300,
        min_segment_length=15000,
        max_segment_length=25000,
        duration=8
    ),
    dna_source=blocks_assembly_station,
    nucleotide_resolution=1000,
    refine_resolution=100,
    a_star_factor=0,
)

numpy.random.seed(1234)
sequence = random_dna_sequence(50000)

import time
t0 = time.time()
quote = chunks_assembly_station.get_quote(sequence,
                                          with_ordering_plan=True)
print time.time() - t0
print (quote)
if quote.accepted:
    print (quote.ordering_plan.summary())

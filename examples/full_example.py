
from dnaadvisor import (ExternalDnaOffer,
                        DnaAssemblyStation,
                        GibsonAssemblyMethod,
                        GoldenGateAssemblyMethod,
                        BuildAGenomeAssemblyMethod,
                        DnaSourcesComparator)
from dnachisel import random_dna_sequence
import dnachisel.constraints as cst
import numpy


# OLIGO COMPANIES

cheap_dna_com = ExternalDnaOffer(
    name="CheapDNA.com",
    sequence_constraints=[lambda seq: "GGTCTC" not in seq,
                          lambda seq: len(seq) < 200],
    price_function=lambda sequence: 0.10 * len(sequence),
    lead_time=10,
    memoize=True
)

deluxe_dna_com = ExternalDnaOffer(
    name="DeluxeDNA.com",
    sequence_constraints=[lambda seq: len(seq) < 200],
    price_function=lambda sequence: 0.20 * len(sequence),
    lead_time=5,
    memoize=True
)

big_dna_com = ExternalDnaOffer(
    name="BigDNA.com",
    sequence_constraints=[lambda seq: 2000 < len(seq) < 4000],
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
    progress_bars=False
)


# BLOCKS TO CHUNKS ASSEMBLY

blocks_sources_comparator = DnaSourcesComparator(
    [
        oligo_assembly_station,
        big_dna_com
    ],
    memoize=True
)


blocks_assembly_comparator = DnaSourcesComparator([
    DnaAssemblyStation(
        name="Blocks Assembly (Gibson)",
        assembly_method=GibsonAssemblyMethod(
            homology_arm_length=40,
            min_segment_length=2000,
            max_segment_length=4000,
            duration=8,
            cost=10
        ),
        dna_source=blocks_sources_comparator,
        nucleotide_resolution=300,
        refine_resolution=False,
        memoize=True,
        progress_bars=True
    ),
    DnaAssemblyStation(
        name="Blocks Assembly (Golden Gate)",
        assembly_method=GoldenGateAssemblyMethod(
            left_overhang='[BsaI]A',
            min_segment_length=2000,
            max_segment_length=4000,
            duration=5,
            cost=2,
            segment_filters=(lambda seq: "GGTCTC" not in seq,)
        ),
        dna_source=blocks_sources_comparator,
        nucleotide_resolution=400,
        refine_resolution=False,
        memoize=True,
        progress_bars=False
    )
])


# CHUNKS TO MEGACHUNKS ASSEMBLY

chunks_assembly_station = DnaAssemblyStation(
    name="Yeast Recombination",
    assembly_method=GibsonAssemblyMethod(
        homology_arm_length=300,
        min_segment_length=15000,
        max_segment_length=25000,
        duration=8
    ),
    dna_source=blocks_assembly_comparator,
    nucleotide_resolution=2000,
    refine_resolution=False,
    progress_bars=True
)


numpy.random.seed(1234)
sequence = random_dna_sequence(50000)

quote = chunks_assembly_station.get_quote(sequence,  with_ordering_plan=True)

print (quote)
if quote.accepted:
    print (quote.ordering_plan.summary())

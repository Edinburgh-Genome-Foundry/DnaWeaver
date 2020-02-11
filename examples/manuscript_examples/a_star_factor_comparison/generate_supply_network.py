import os
from dnaweaver.biotools import gc_content
import dnaweaver as dw


def generate_supply_network(a_star_factor):

    oligo_com = dw.CommercialDnaOffer(
        name="Oligo.com",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=200)],
        pricing=dw.PerBasepairPricing(0.10),
        lead_time=7,
    )

    deluxe_dna_com = dw.CommercialDnaOffer(
        name="DeluxeDNA.com",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=4000)],
        pricing=dw.PerBasepairPricing(0.20),
        lead_time=10,
    )

    cheap_dna_com = dw.CommercialDnaOffer(
        name="CheapDNA.com",
        sequence_constraints=[
            dw.SequenceLengthConstraint(max_length=4000),
            dw.NoPatternConstraint(enzyme="AarI"),
            dw.NoPatternConstraint(enzyme="BsaI"),
            lambda seq: (0.4 < gc_content(seq) < 0.6),
        ],
        pricing=dw.PerBasepairPricing(0.10),
        lead_time=15,
    )

    # OLIGOS TO BLOCKS ASSEMBLY

    oligo_assembly_station = dw.DnaAssemblyStation(
        name="Oligo Assembly Station",
        assembly_method=dw.OligoAssemblyMethod(
            overhang_selector=dw.TmSegmentSelector(
                min_size=15, max_size=25, min_tm=50, max_tm=70
            ),
            min_segment_length=40,
            max_segment_length=200,
            sequence_constraints=[dw.SequenceLengthConstraint(max_length=1500)],
            duration=8,
            cost=2,
        ),
        supplier=oligo_com,
        coarse_grain=20,
        fine_grain=False,
        a_star_factor=a_star_factor,
    )

    # BLOCKS TO CHUNKS ASSEMBLY

    blocks_sources_comparator = dw.DnaSuppliersComparator(
        name="bs_comparator",
        suppliers=[oligo_assembly_station, cheap_dna_com, deluxe_dna_com],
        memoize=True,
    )

    gibson_blocks_assembly_station = dw.DnaAssemblyStation(
        name="Gibson Blocks Assembly",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.FixedSizeSegmentSelector(80),
            min_segment_length=1000,
            max_segment_length=4000,
            duration=8,
            cost=16,
        ),
        supplier=blocks_sources_comparator,
        coarse_grain=300,
        fine_grain=False,
        memoize=True,
        a_star_factor=a_star_factor,
    )

    goldengate_blocks_assembly_station = dw.DnaAssemblyStation(
        name="Golden Gate Blocks Assembly",
        assembly_method=dw.GoldenGateAssemblyMethod(
            enzyme="BsmBI",
            wildcard_basepair="A",
            min_segment_length=1000,
            max_segment_length=4000,
            duration=5,
            cost=6,
        ),
        supplier=blocks_sources_comparator,
        coarse_grain=400,
        fine_grain=False,
        memoize=True,
        a_star_factor=a_star_factor,
    )
    ecoli_genome_path = os.path.join(
        "..", "..", "data", "ecoli_blast_db", "ecoli"
    )
    ecoli_genome = dw.PcrExtractionStation(
        "E. coli Genome (PCR)",
        primers_supplier=oligo_com,
        homology_selector=dw.TmSegmentSelector(),
        blast_database=ecoli_genome_path,
        max_amplicon_length=10000,
        extra_time=3,
        extra_cost=1,
    )

    # CHUNKS TO MEGACHUNKS ASSEMBLY

    return dw.DnaAssemblyStation(
        name="Chunks assembly (Yeast)",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.FixedSizeSegmentSelector(300),
            min_segment_length=7000,
            max_segment_length=25000,
            duration=8,
        ),
        supplier=[
            ecoli_genome,
            goldengate_blocks_assembly_station,
            gibson_blocks_assembly_station,
        ],
        coarse_grain=1000,
        fine_grain=None,
        a_star_factor=a_star_factor,
        memoize=True,
    )

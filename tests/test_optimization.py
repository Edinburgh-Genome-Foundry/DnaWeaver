import os
from dnaweaver import (
    CommercialDnaOffer,
    GibsonAssemblyMethod,
    PerBasepairPricing,
    DnaAssemblyStation,
    NoPatternConstraint,
    FixedSizeSegmentSelector,
    DnaSuppliersComparator,
    NoPatternConstraint,
    SequenceLengthConstraint,
)
from dnachisel import (
    DnaOptimizationProblem,
    load_record,
    EnforceGCContent,
    EnforceTranslation,
)
from dnaweaver.utils import OptimizeManufacturability


def test_optimization_1():
    company_ingen = CommercialDnaOffer(
        name="Company InGen",
        pricing=PerBasepairPricing(0.08),
        sequence_constraints=[NoPatternConstraint(enzyme="AarI")],
    )
    company_delux = CommercialDnaOffer(
        name="Company Delux",
        pricing=PerBasepairPricing(0.66),
        sequence_constraints=[],
    )

    assembly_station = DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=GibsonAssemblyMethod(
            overhang_selector=FixedSizeSegmentSelector(20),
            min_segment_length=200,
            max_segment_length=1200,
        ),
        supplier=[company_ingen, company_delux],
        coarse_grain=20,
        # a_star_factor="auto",
    )
    sequence_path = os.path.join(
        "tests", "data", "test_optimization_sequence_1.fa"
    )
    sequence = load_record(sequence_path)
    objective = OptimizeManufacturability(assembly_station)
    problem = DnaOptimizationProblem(sequence=sequence, objectives=[objective])
    quote = objective.get_quote(problem)
    score = problem.objective_scores_sum()
    assert -367 < score < -366
    problem.randomization_threshold = 0
    problem.max_random_iters = 5
    problem.optimize()
    score = problem.objective_scores_sum()
    assert -244 < score < -243


def test_optimization_2():
    sequence_path = os.path.join(
        "tests", "data", "test_optimization_sequence_2.fa"
    )
    sequence = str(load_record(sequence_path).seq)[:5500]

    deluxe_dna = CommercialDnaOffer(
        name="DeluxeDNA.com",
        sequence_constraints=[SequenceLengthConstraint(max_length=4000)],
        pricing=PerBasepairPricing(0.20),
        lead_time=10,
    )

    cheap_dna = CommercialDnaOffer(
        name="CheapDNA.com",
        sequence_constraints=[
            NoPatternConstraint(enzyme="BsaI"),
            EnforceGCContent(0.3, 0.7, window=60),
        ],
        pricing=PerBasepairPricing(0.10),
        lead_time=15,
    )

    # BLOCKS TO CHUNKS ASSEMBLY

    gibson_blocks_assembly_station = DnaAssemblyStation(
        name="Gibson Blocks Assembly",
        assembly_method=GibsonAssemblyMethod(
            overhang_selector=FixedSizeSegmentSelector(10),
            min_segment_length=1000,
            max_segment_length=6000,
            duration=8,
            cost=16,
        ),
        supplier=[deluxe_dna, cheap_dna],
        coarse_grain=30,
        fine_grain=False,
        memoize=True,
        # a_star_factor="auto",
    )

    quote_before = gibson_blocks_assembly_station.get_quote(sequence)
    assert quote_before.price > 850

    objective = OptimizeManufacturability(gibson_blocks_assembly_station)

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation(location=(0, 4998))],
        objectives=[objective],
    )

    problem.randomization_threshold = 0  # Forces "random search" mode
    problem.max_random_iters = 5
    problem.optimize()

    print("OPTIMIZATION DONE, GENERATING REPORT")

    quote_after = gibson_blocks_assembly_station.get_quote(problem.sequence)
    assert quote_after.price < 580

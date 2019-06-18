import os
from dnaweaver import (CommercialDnaOffer, GibsonAssemblyMethod,
                       PerBasepairPricing, DnaAssemblyStation,
                       NoPatternConstraint, FixedSizeOverhangSelector)
from dnachisel import DnaOptimizationProblem, load_record
from dnaweaver.OptimizeManufacturability import OptimizeManufacturability

def test_optimization_1():
    company_ingen = CommercialDnaOffer(
        name="Company InGen",
        pricing=PerBasepairPricing(0.08),
        sequence_constraints=[NoPatternConstraint(enzyme='AarI')]
    )
    company_delux = CommercialDnaOffer(
        name="Company Delux",
        pricing=PerBasepairPricing(0.66),
        sequence_constraints=[],
    )

    assembly_station = DnaAssemblyStation(
        name='Gibson Assembly Station',
        assembly_method=GibsonAssemblyMethod(
            overhang_selector=FixedSizeOverhangSelector(20),
            min_segment_length=200,
            max_segment_length=1200,
        ),
        supplier=[company_ingen, company_delux],
        coarse_grain=20,
        a_star_factor = 'auto'
    )
    sequence_path = os.path.join('tests', 'test_optimization_sequence.fa')
    sequence = load_record(sequence_path)
    objective = OptimizeManufacturability(assembly_station)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        objectives=[objective]
    )
    quote = objective.get_quote(problem)
    score = problem.objective_scores_sum()
    assert -367 < score < -366
    problem.randomization_threshold = 0
    problem.max_random_iters = 5
    problem.optimize()
    quote = objective.get_quote(problem)
    score = problem.objective_scores_sum()
    assert -244 < score < -243
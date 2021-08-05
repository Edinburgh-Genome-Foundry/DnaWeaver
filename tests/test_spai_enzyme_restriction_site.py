"""Test contributed by jlerman44/alexsongrj to add SapI enzyme restriction sites #10"""
import dnaweaver as dw
from dnaweaver import TmSegmentSelector
import dnacauldron as dc
from dnachisel import (
    random_dna_sequence,
    DnaOptimizationProblem,
    AvoidPattern,
    AvoidChanges,
)
import pandas as pd
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from geneblocks import sequences_are_circularly_equal

restriction_enzyme = "SapI"


def test_golden_gate_assembly(tmpdir):

    desired_sequence = ""
    sequences = {
        "fragment 1": random_dna_sequence(1000, seed=123),
        "fragment 2": random_dna_sequence(1000, seed=456),
        "fragment 3": random_dna_sequence(1000, seed=789),
    }
    for fragment in sequences:
        sequence = sequences[fragment]
        desired_sequence += sequence
    # make sure desired has no SapI site.
    constraints = []
    constraints.append(AvoidPattern("{}_site".format(restriction_enzyme)))
    problem = DnaOptimizationProblem(
        sequence=desired_sequence,
        constraints=constraints,
        objectives=[AvoidChanges().as_passive_objective()],
    )
    problem.resolve_constraints()
    problem.optimize()
    domesticated_sequence = problem.sequence
    desired_sequence = dw.SequenceString(
        domesticated_sequence, metadata={"topology": "circular"}
    )

    oligos_company = dw.CommercialDnaOffer(
        "OligoCompany",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=200)],
        pricing=dw.PerBasepairPricing(0.1),
    )
    pcr_station = dw.PcrExtractionStation(
        name="PCR station",
        max_overhang_length=50,
        primers_supplier=oligos_company,
        sequences=sequences,
        homology_selector=TmSegmentSelector(min_tm=35, max_tm=70),
        extra_cost=5,
    )
    assembly_station = dw.DnaAssemblyStation(
        name="Golden Gate assembly",
        assembly_method=dw.GoldenGateAssemblyMethod(
            enzyme=restriction_enzyme,
            min_segment_length=100,
            max_segment_length=10000,
            wildcard_basepair="G",
            left_addition="TCCTAG",
            right_addition="CTAGGA",
        ),
        supplier=pcr_station,
        coarse_grain=10,
        fine_grain=None,
        logger="bar",
    )

    # THIS LINE WILL PRE-BLAST THE SEQUENCE TO ACCELERATE COMPUTATIONS.
    assembly_station.prepare_network_on_sequence(desired_sequence)

    # FIND AN ASSEMBLY PLAN AND PRINT IT.
    quote = assembly_station.get_quote(desired_sequence)
    assembly_report = quote.to_assembly_plan_report()
    target = os.path.join(str(tmpdir), "solution")
    os.mkdir(target)
    assert os.listdir(target) == []
    assembly_report.write_full_report(target)
    assert os.listdir(target) != []

    # now verify things are correct
    df = pd.read_csv(f"{target}/sequences.csv", sep=";")
    target_sequence = df[df.Source == "Golden Gate assembly"].Sequence.unique()[0]
    df = df[df.Source == "PCR station"]

    parts = {}
    ids = []
    for _, row in df.iterrows():
        new_seq_record = SeqRecord(Seq(row.Sequence))
        new_seq_record.id = row.ID
        new_seq_record.description = row.ID
        new_seq_record.name = row.ID
        new_seq_record.annotations["topology"] = "linear"
        parts[row.ID] = new_seq_record
        ids.append(row.ID)

    repository = dc.SequenceRepository(collections={"my_collection": parts})

    assembly = dc.Type2sRestrictionAssembly(parts=ids, enzyme=restriction_enzyme)
    simulation = assembly.simulate(sequence_repository=repository)

    assert len(simulation.construct_records) == 1
    predicted_construct_record = simulation.construct_records[0]

    assert sequences_are_circularly_equal([predicted_construct_record, target_sequence])

import os
import dnaweaver as dw
import dnacauldron as dc
from geneblocks import sequences_are_circularly_equal

this_directory = os.path.dirname(os.path.realpath(__file__))
plasmid_path = os.path.join(this_directory, "circular_no_BsmBI.gb")

def extract_records_from_quote(quote):
    quote.compute_full_assembly_plan()
    records = {
        q.id: dc.sequence_to_biopython_record(q.sequence, id=q.id)
        for loc, q in quote.assembly_plan.items()
    }
    repo = dc.SequenceRepository(collections={'parts': records})
    record_names = list(records.keys())
    return record_names, repo

def test_circular_golden_gate():

    record = dw.load_record(plasmid_path)
    sequence = dw.SequenceString.from_record(record)
    assert sequence.metadata['topology'] == "circular"
    assert dw.get_sequence_topology(sequence) == "circular"

    commercial_provider = dw.CommercialDnaOffer(
        name="supplier", pricing=dw.PerBasepairPricing(0.1),
    )
    assembly_station = dw.DnaAssemblyStation(
        name="golden_gate",
        supplier=commercial_provider,
        assembly_method=dw.GoldenGateAssemblyMethod(
            enzyme="BsmBI",
            max_segment_length=3000
        ),
        coarse_grain=600,
        fine_grain=None,
        cut_spread_radius=2,
    )
    quote = assembly_station.get_quote(sequence)
    assert quote.accepted
    assert 1200 < quote.price < 1220

    # ROUND TRIP: SIMULATE CLONING TO CHECK ASSEMBLY PLAN VALIDITY

    part_names, repo = extract_records_from_quote(quote)
    quote = assembly_station.get_quote(sequence)
    assembly = dc.Type2sRestrictionAssembly(parts=part_names)
    simulation = assembly.simulate(repo)
    assert len(simulation.construct_records) == 1
    simulated_record = simulation.construct_records[0]
    assert sequences_are_circularly_equal([record, simulated_record])

def test_circular_gibson():

    record = dw.load_record(plasmid_path)
    sequence = dw.SequenceString.from_record(record)
    assert sequence.metadata['topology'] == "circular"
    assert dw.get_sequence_topology(sequence) == "circular"

    commercial_provider = dw.CommercialDnaOffer(
        name="supplier", pricing=dw.PerBasepairPricing(0.1),
    )
    assembly_station = dw.DnaAssemblyStation(
        name="golden_gate",
        supplier=commercial_provider,
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.TmSegmentSelector(),
            max_segment_length=3000
        ),
        coarse_grain=600,
        fine_grain=None,
        cut_spread_radius=2,
    )
    quote = assembly_station.get_quote(sequence)
    assert quote.accepted
    assert 1200 < quote.price < 1220

    # ROUND TRIP: SIMULATE CLONING TO CHECK ASSEMBLY PLAN VALIDITY

    part_names, repo = extract_records_from_quote(quote)
    quote = assembly_station.get_quote(sequence)
    assembly = dc.GibsonAssembly(parts=part_names)
    simulation = assembly.simulate(repo)
    assert len(simulation.construct_records) == 1
    simulated_record = simulation.construct_records[0]
    print (len(record), len(simulated_record))
    assert sequences_are_circularly_equal([record, simulated_record])

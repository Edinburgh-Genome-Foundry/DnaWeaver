# This example feeds 2 sequences, one with BsmBI sites, the other without, to
# a comparator Gibson/Golden Gate assembly.
# For more realism, the folling adapters as added to the network:
# - A Golden Gate adapter which adds flanks TAGG-ACGA to the sequence for
#   assembly on an EMMA plasmid
# - A Gibson Assembly adapter which adds flanks to the sequence for assembly on
#   a backbone
# - A linearization station which linearizes parts provided on a backbone by
#   the commercial provider Plasmideroo.

import dnaweaver as dw
import os

this_directory = os.path.dirname(os.path.realpath(__file__))


def test_example_with_adapters():
    oligo_com = dw.CommercialDnaOffer(
        name="Oligo.com",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=200)],
        pricing=dw.PerBasepairPricing(0.10),
        lead_time=7,
    )

    plasmideroo = dw.CommercialDnaOffer(
        name="Plasmideroo",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=4000)],
        pricing=dw.PerBasepairPricing(0.20),
        lead_time=10,
    )

    plasmideroo_linearization = dw.PcrLinearizationStation(
        name="Gibson Linearization",
        supplier=plasmideroo,
        primers_supplier=oligo_com,
        homology_selector=dw.TmSegmentSelector(),
    )

    gibson_blocks_assembly_station = dw.DnaAssemblyStation(
        name="Gibson Assembly",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.FixedSizeSegmentSelector(80),
            min_segment_length=1000,
            max_segment_length=4000,
            duration=8,
            cost=16,
        ),
        supplier=plasmideroo_linearization,
        coarse_grain=100,
        fine_grain=False,
        a_star_factor="auto",
    )

    gibson_station_adapter = dw.SequenceAdapter(
        name="Gibson Adapter",
        supplier=gibson_blocks_assembly_station,
        left_addition="TAAGACTAGTGCACATAATACGG",
        right_addition="ATTGTCACACTACAAACATGAC",
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
        supplier=plasmideroo,
        coarse_grain=100,
        fine_grain=False,
        a_star_factor="auto",
    )

    golden_gate_station_adapter = dw.SequenceAdapter(
        name="Golden Gate Adapter",
        supplier=goldengate_blocks_assembly_station,
        left_addition="TAGG",
        right_addition="ACGA",
    )

    main_source = dw.DnaSuppliersComparator(
        [gibson_station_adapter, golden_gate_station_adapter]
    )

    # RESOLUTION

    print("\n\nSEQUENCE WITHOUT BMSBI SITES:\n")
    path = os.path.join(this_directory, "sequence_no_bsmbi.fa")
    sequence = dw.load_record(path)
    quote = main_source.get_quote(sequence, with_assembly_plan=True)
    assert quote.source == goldengate_blocks_assembly_station
    print(quote.assembly_step_summary())
    assert 2000 < quote.price < 2020

    print("\n\nSEQUENCE WITH BMSBI SITES:\n")
    path = os.path.join(this_directory, "sequence_bsmbi.fa")
    sequence = dw.load_record(path)
    quote = main_source.get_quote(sequence, with_assembly_plan=True)
    print(quote.assembly_step_summary())
    assert quote.source == gibson_blocks_assembly_station
    assert 2050 < quote.price < 2100
    
    # FOR COVERAGE

    assembly_report = quote.to_assembly_plan_report()
    assembly_report.write_full_report('@memory')

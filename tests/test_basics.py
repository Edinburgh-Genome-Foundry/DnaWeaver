from dnaweaver import (
    DnaAssemblyStation,
    GibsonAssemblyMethod,
    TmSegmentSelector,
    CommercialDnaOffer,
    PerBasepairPricing,
    random_dna_sequence,
)


def test_simple_gibson_assembly_station():
    dna_provider = CommercialDnaOffer(
        name="Company InGen", pricing=PerBasepairPricing(0.08)
    )
    assembly_station = DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=GibsonAssemblyMethod(
            overhang_selector=TmSegmentSelector(),
            min_segment_length=300,
            max_segment_length=1200,
        ),
        supplier=dna_provider,
        coarse_grain=10,
    )
    sequence = random_dna_sequence(5000, seed=1234)
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
    assert quote.accepted
    assert 405.7 < quote.price < 405.8

    # Now with a lower-case sequence
    quote = assembly_station.get_quote(sequence.lower(), with_assembly_plan=True)
    assert quote.accepted
    assert 405.7 < quote.price < 405.8

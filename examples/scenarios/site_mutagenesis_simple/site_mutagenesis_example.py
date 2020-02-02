import dnaweaver as dw
import os

oligos_company = dw.CommercialDnaOffer(
    "OligoCompany",
    sequence_constraints=[dw.SequenceLengthConstraint(max_length=200)],
    pricing=dw.PerBasepairPricing(0.1)
)
ecoli_db_path = os.path.join('..', '..', 'data', 'ecoli_blast_db', 'ecoli')
pcr_station = dw.PcrExtractionStation(
    name="PCR station",
    max_overhang_length=50,
    homology_selector=dw.TmSegmentSelector(),
    primers_supplier=oligos_company,
    blast_database=ecoli_db_path,
    extra_cost=5
)
assembly_station = dw.DnaAssemblyStation(
    name="Golden Gate assembly",
    assembly_method=dw.GoldenGateAssemblyMethod(enzyme='BsaI'),
    supplier=pcr_station,
    coarse_grain=100,
    fine_grain=0,
    logger='bar'
)

# LOAD THE (SITE-FREE) DESIRED SEQUENCE
desired_sequence = str(dw.load_record("desired_sequence.gb").seq)

# THIS LINE WILL PRE-BLAST THE SEQUENCE TO ACCELERATE COMPUTATIONS.
assembly_station.prepare_network_on_sequence(desired_sequence)

# FIND AN ASSEMBLY PLAN AND PRINT IT.
quote = assembly_station.get_quote(desired_sequence, with_assembly_plan=True)
print(quote.assembly_step_summary())

quote.compute_fragments_final_locations()
report = quote.to_assembly_plan_report()
report.write_full_report("report.zip")
print ("Done! (see report.zip)")

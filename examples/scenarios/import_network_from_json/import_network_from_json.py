import os
import dnaweaver as dw
import json
ecoli_db_path = os.path.join("..", "..", "data", "ecoli_blast_db", "ecoli")
dw.PcrExtractionStation.dna_banks = {"e_coli": ecoli_db_path}
dw.PartsLibrary.collections_by_id["EMMA"] = dict(
    library_class="golden_gate", fasta_file="emma_parts.fa"
)
sequence = dw.load_record("sequence_with_emma_and_ecoli_parts.fa")
main_station = dw.DnaSupplier.from_json_data("state_from_web_app.json")
main_station.prepare_network_on_sequence(sequence)
quote = main_station.get_quote(sequence)
print(quote.assembly_step_summary())

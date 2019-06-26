from dnaweaver import (supply_network_from_json,
                       PartsLibrary,
                       PcrExtractionStation,
                       load_record)
import json
import os
with open("sequence_with_emma_and_ecoli_parts.txt", "r") as f:
    sequence = f.read()
with open('state_from_web_app.json', 'r') as f:
    data = json.load(f)

PcrExtractionStation.dna_banks = {
    'e_coli': os.path.join("..", "..", "data",
                           "ecoli_blast_db", "ecoli"),
}
PartsLibrary.collections_by_id['EMMA'] = dict(
    library_class='golden_gate',
    fasta_file="emma_parts.fa"
)
levels, network, main_id = supply_network_from_json(data['graph'])
main = network[main_id]
main.prepare_network_on_sequence(sequence)
quote = main.get_quote(sequence)
print (quote.assembly_step_summary())

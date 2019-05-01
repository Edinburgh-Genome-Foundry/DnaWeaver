from base64 import b64decode
import json
import os
from dnaweaver import supply_network_from_json, PcrOutStation

THIS_DIR = os.path.join("tests", "test_scenarios",
                        "test_domestication_from_json")
PcrOutStation.dna_banks = {
    'e_coli': os.path.join('tests', 'data', 'ecoli_blast_db', 'ecoli')
}

def test_supply_network_from_json():
    with open(os.path.join(THIS_DIR, "domestication.json"), "r") as f:
        data = json.load(f)
    levels, suppliers_dict, main_id = supply_network_from_json(data['graph'])
    main = suppliers_dict[main_id]
    sequence_content = data['sequence_file']['content']
    sequence = b64decode(sequence_content.split('base64,')[1]).decode()
    main.prepare_network_on_sequence(sequence)
    quote = main.get_quote(sequence)
    assert 130 < quote.price < 139
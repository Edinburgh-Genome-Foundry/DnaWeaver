from base64 import b64decode
import json
import os
from dnaweaver import DnaSupplier, PcrExtractionStation

THIS_DIR = os.path.join(
    "tests", "test_scenarios", "test_domestication_from_json"
)
PcrExtractionStation.dna_banks = {
    "e_coli": os.path.join("tests", "data", "ecoli_blast_db", "ecoli")
}


def test_supply_network_from_json():
    path = os.path.join(THIS_DIR, "domestication.json")
    with open(path, "r") as f:
        data = json.load(f)
    main = DnaSupplier.from_json_data(data=data['graph'])
    sequence_content = data["sequence_file"]["content"]
    sequence = b64decode(sequence_content.split("base64,")[1]).decode()
    main.prepare_network_on_sequence(sequence)
    quote = main.get_quote(sequence)
    assert 130 < quote.price < 160

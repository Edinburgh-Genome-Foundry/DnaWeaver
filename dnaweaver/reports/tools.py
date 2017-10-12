from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


class ObjectDict(dict):

    def __getattr__(self, key):
        return self.__dict__[key]

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    @staticmethod
    def from_dict(d):
        obj = ObjectDict({
            key: ([
                ObjectDict.from_dict(e)
                if isinstance(e, dict) else e
                for e in value
            ]
                if isinstance(value, (list, tuple))
                else (ObjectDict.from_dict(value)
                      if isinstance(value, dict)
                      else value))
            for (key, value) in d.items()
        })
        for key, value in obj.items():
            obj.__dict__[key] = value
        return obj

class JsonQuote:

    def __init__(self, tree, sources):
        self.tree = ObjectDict.from_dict(tree)
        self.sources = sources

    def to_quotes_list(self):
        tree = deepcopy(self.tree)
        nodes = []

        def rec(node, depth=0):
            if node.get("_visited", False):
                return
            node["_visited"] = True
            assembly_plan = node.pop("assembly_plan")
            node["children"] = [
                n["id"] for n in node.get("assembly_plan", [])
            ]
            nodes.append(node)
            for other in sorted(assembly_plan,
                                key=lambda n: n["segment_start"]):
                rec(other)
        rec(tree)
        return nodes

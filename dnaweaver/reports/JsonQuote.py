from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

try:  # Python2 vs. Python3 ways of dealing with file-like strings
    from StringIO import StringIO
    USE_BYTES = False
except ImportError:
    from io import StringIO, BytesIO
    USE_BYTES = True

class ObjectDict(dict):

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        dict.__getattribute__(self, key)

    def __setattr__(self, key, value):
        self[key] = value
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
            sanitized_key = key.replace(" ", "_").replace(".", "_")
            obj.__dict__[sanitized_key] = value
        return obj


class JsonQuote:

    def __init__(self, tree, sources):
        self.tree = ObjectDict.from_dict(tree)
        self.sources = ObjectDict.from_dict(sources)

    @staticmethod
    def from_dnaweaver_quote(quote):
        tree = quote.assembly_tree_as_dict()
        sources = quote.source.dict_supply_graph()
        return JsonQuote(tree, sources)

    def to_quotes_list(self):
        tree = deepcopy(self.tree)
        nodes = []

        def rec(node, depth=0):
            if node.get("_visited", False):
                return
            node["_visited"] = True
            assembly_plan = node.get("assembly_plan", [])
            node["children"] = [
                n["id"] for n in assembly_plan
            ]
            nodes.append(node)
            for other in sorted(assembly_plan,
                                key=lambda n: n["segment_start"]):
                rec(other)
        rec(tree)
        return nodes

    def quote_tree_to_SeqRecord(self, record=None, record_id=None):
        """Return a Biopython seqrecord of the quote.

        >>> record = to_SeqRecord(solution)
        >>> # Let's plot with DnaVu:
        >>> from dnavu import create_record_plot
        >>> from bokeh.io import output_file, show
        >>> output_file("view.html")
        >>> plot = create_record_plot(record)
        >>> show(plot)
        """
        if record_id is None:
            record_id = self.id
        if record is None:
            record = SeqRecord(Seq(self.sequence, DNAAlphabet()), id=record_id)
        else:
            record = deepcopy(record)

        if self.assembly_plan is not None:
            features = [
                SeqFeature(
                    FeatureLocation(segment[0], segment[1], 1),
                    type="Feature",
                    qualifiers={
                        "label": "%s - From %s" % (quote.id, quote.source),
                        "name": quote.id,
                        "source": quote.source,
                        "price": quote.price,
                        "lead_time": quote.lead_time
                    }
                )
                for segment, quote in self.assembly_plan.items()
            ]
            record.features = features + record.features
        return record

    def quote_tree_to_genbank(self, filename=None, filehandle=None,
                              record=None, record_id=None):
        record = self.to_SeqRecord(record=record, record_id=record_id)
        if filename is not None:
            with open(filename, "w+") as f:
                SeqIO.write(record, f, "genbank")
        else:
            output = StringIO()
            SeqIO.write(record, output, "genbank")
            return output.getvalue()

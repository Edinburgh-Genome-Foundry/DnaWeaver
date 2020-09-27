import json
from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from io import StringIO

from ..AssemblyPlanReport import AssemblyPlanReport


class ExportsMixin:
    def tree_as_list(self):
        """Return a list containing the current AssemblyOperation and all its
        sub-operations and their respective sub-operations.

        Said otherwise, it flattens the assembly tree into the list of all
        nodes.
        """
        result = [self]
        if self.assembly_plan is not None:
            result += sum(
                [child.tree_as_list() for segment, child in self.assembly_plan.items()],
                [],
            )
        return result

    def assembly_plan_as_dict(self, as_json=False, json_indent=None):
        """Return a JSON-like version of the nested tree.

        Parameters
        ----------

        as_json
          If True, a JSON string is returned, else the result is a dict object.

        json_indent
          number of spaces in the JSON indentation (for pretty printing). The
          default None means that the JSON will be on one line (TODO: check).


        Returns
        -------

          {
          "id": self.id,
          "source": self.source.name,
          "price": self.price,
          "lead_time": self.lead_time,
          "sequence": self.sequence,
          "message": self.message,
          "metadata" = self.metadata,
          "assembly_plan": { (start1, end1): {(subquote_1)},
                             (start2, end2): {(subquote_2)},
                           }
          }
        """
        final_location = (
            self.final_location if hasattr(self, "final_location") else None
        )
        matching_segment = (
            self.matching_segment if hasattr(self, "matching_segment") else None
        )

        assembly_plan = []
        if self.assembly_plan is not None:
            for (segment, quote) in self.assembly_plan.items():
                quote_as_dict = quote.assembly_plan_as_dict()
                quote_as_dict["segment_start"] = segment[0]
                quote_as_dict["segment_end"] = segment[1]
                assembly_plan.append(quote_as_dict)
        tree = {
            "id": self.id,
            "source": self.source.name,
            "price": self.price,
            "lead_time": self.lead_time,
            "sequence": self.sequence,
            "message": self.message,
            "metadata": self.metadata,
            "assembly_plan": assembly_plan,
            "final_location": final_location,
            "matching_segment": matching_segment,
            "accepted": self.accepted,
        }
        metadata = tree["metadata"]
        if "via" in metadata:
            metadata["via"] = [
                station if isinstance(station, str) else station.name
                for station in metadata["via"]
            ]

        if as_json:
            return json.dumps(tree, indent=json_indent)
        else:
            return tree

    def to_record(self, record=None, record_id=None):
        """Return a Biopython seqrecord of the quote.

        >>> record = to_record(solution)
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
            if has_dna_alphabet:  # Biopython <1.78
                record = SeqRecord(Seq(self.sequence, DNAAlphabet()), id=record_id)
            else:
                record = SeqRecord(Seq(self.sequence), id=record_id)
            record.annotations["molecule_type"] = "DNA"
        else:
            record = deepcopy(record)

        if self.assembly_plan is not None:
            features = [
                SeqFeature(
                    FeatureLocation(segment[0], segment[1], 1),
                    type="Feature",
                    qualifiers={
                        "name": quote.id,
                        "source": quote.source,
                        "price": quote.price,
                        "lead_time": quote.lead_time,
                    },
                )
                for segment, quote in self.assembly_plan.items()
            ]
            record.features = features + record.features
        return record

    def write_genbank(
        self, filename=None, filehandle=None, record=None, record_id=None
    ):
        record = self.to_record(record=record, record_id=record_id)
        if filename is not None:
            with open(filename, "w+") as f:
                SeqIO.write(record, f, "genbank")
        else:
            output = StringIO()
            SeqIO.write(record, output, "genbank")
            return output.getvalue()

    def to_assembly_plan_report(
        self, refine_fragments_locations=True, autocolor_quotes=True
    ):
        """Convert the quote into a full assembly plan data structure which
        can be used to generate assembly reports."""
        if refine_fragments_locations:
            self.compute_fragments_final_locations()
        if not self.full_assembly_plan_computed:
            self.compute_full_assembly_plan()
        original_source = self.source
        if "via" in self.metadata:
            # intermediary comparator of the quote
            original_source = self.metadata["via"][0]
        report = AssemblyPlanReport(
            plan=self.assembly_plan_as_dict(),
            sources=original_source.dict_supply_graph(),
        )
        if autocolor_quotes:
            report.autocolor_quote_sources()
        return report

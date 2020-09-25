from copy import deepcopy
from io import StringIO
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


class GenbankExportMixin:
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

        if self.plan is not None:
            features = [
                SeqFeature(
                    FeatureLocation(q.segment_start, q.segment_end, 1),
                    type="misc_feature",
                    qualifiers={
                        "label": "%s - From %s" % (q.id, q.source),
                        "name": q.id,
                        "source": q.source,
                        "price": q.price,
                        "lead_time": q.lead_time,
                    },
                )
                for q in self.plan
            ]
            record.features = features + record.features
        return record

    def write_genbank(
        self, filename=None, filehandle=None, record=None, record_id=None
    ):
        record = self.to_record(record=record, record_id=record_id)
        if filehandle is None:
            with open(filename, "w+") as f:
                SeqIO.write(record, f, "genbank")
        else:
            SeqIO.write(record, filehandle, "genbank")

    def write_all_sequence_records(self, target):
        for step in self.to_steps_list():
            record = self.plan_step_to_record(step)
            path = target._file(step.id + ".gb").open("w")
            SeqIO.write(record, path, "genbank")

    @staticmethod
    def plan_step_to_record(plan_step, record=None, record_id=None):
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
            record_id = plan_step.id
        if record is None:
            if has_dna_alphabet:  # Biopython <1.78
                record = SeqRecord(Seq(plan_step.sequence, DNAAlphabet()), id=record_id)
            else:
                record = SeqRecord(Seq(plan_step.sequence), id=record_id)
            record.annotations["molecule_type"] = "DNA"
        else:
            record = deepcopy(record)

        if plan_step.assembly_plan is not None:
            features = [
                SeqFeature(
                    FeatureLocation(q.segment_start, q.segment_end, 1),
                    type="misc_feature",
                    qualifiers={
                        "label": "%s - From %s" % (q.id, q.source),
                        "name": q.id,
                        "source": q.source,
                        "price": q.price,
                        "lead_time": q.lead_time,
                    },
                )
                for q in plan_step.assembly_plan
            ]
            record.features = features + record.features
        return record

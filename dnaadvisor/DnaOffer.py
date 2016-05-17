from copy import deepcopy
from dnachisel.constraints import Constraint
from dnachisel import DnaCanvas

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except:
    PANDAS_AVAILABLE = False

try:
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import DNAAlphabet
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    BIOPYTHON_AVAILABLE = True
except:
    BIOPYTHON_AVAILABLE = False


class DnaOfferEvaluation:
    """Evaluates an offer on a sequence.

    Parameter
    ---------

    sequence
      An ATGC string of a DNA sequence

    offer
      a DnaOffer object (with constraints and pricing)

    segment
      Optional. A couple (start, end) defining the segment between two "cuts"
      in the original sequence on which that segment was considered.

    """

    def __init__(self, sequence, offer, segment=None):

        self.sequence = sequence
        self.offer = offer
        self.segment = segment

        # The next paragraph is there
        constraints = offer.constraints
        dnachisel_constraints = [
            constraint for constraint in constraints
            if isinstance(constraint, Constraint)
        ]

        if dnachisel_constraints != []:
            canvas = DnaCanvas(sequence, dnachisel_constraints)
            constraints = [
                constraint for constraint in constraints
                if not isinstance(constraint, Constraint)
            ] + [lambda seq: canvas.all_constraints_pass()]
        if not all(constraint(sequence) for constraint in constraints):
            self.is_orderable = False
            self.price = None
        else:
            self.is_orderable = True
            self.price = offer.pricing(sequence)
            self.basepair_price = 1.0* self.price / len(self.sequence)

    def __repr__(self):
        if not self.is_orderable:
            return "%s (not orderable)" % self.offer.name
        else:
            return "%s%s %.02f$" % ("" if self.segment is None
                                    else ("%s " % str(self.segment)),
                                    self.offer.name, self.price)


class DnaOffer:

    def __init__(self, name, constraints, pricing):
        self.name = name
        self.constraints = constraints
        self.pricing = pricing

    def __repr__(self):
        return self.name


class DnaOrderingPlan:
    """Class for representing a group of orders and exporting to many formats.

    Parameters
    ----------

    plan
      A dictionary { sequence: offer_evaluation, ... } where `sequence`
      is an ATGC string sequence of DNA, and offer_evaluation is a
      DnaOfferEvaluation.

    """

    def __init__(self, offer_evaluations, full_sequence=None):
        self.plan = offer_evaluations
        self.full_sequence = full_sequence

    def summary(self):
        """Return a print-friendly, simple string of the ordering plan."""
        evaluations = sorted(self.plan, key=lambda o: o.segment)
        plan = "\n  ".join(str(e) for e in evaluations)
        return "Ordering plan:\n  %s\n  Total:%d$" % (plan, self.total_price())

    def total_price(self):
        return sum(evaluation.price for evaluation in self.plan)

    def cuts_indices(self):
        return sorted([offer.segment[1] for offer in self.plan])[:-1]

    def write_pdf_report(self):
        """Write a PDF report of the ordering plan."""
        raise NotImplementedError()

    def to_html_widget(self):
        """Return an interactive widget of the ordering"""
        raise NotImplementedError()

    def to_dataframe(self, sort_by="segment"):
        """Export as a Pandas dataframe. Particularly useful for Notebooks"""

        if not PANDAS_AVAILABLE:
            raise ImportError("Install Pandas to use to_dataframe.")
        dataframe = pd.DataFrame.from_records([
            {"offer": evaluation.offer,
             "segment": evaluation.segment,
             "length": evaluation.segment[1] - evaluation.segment[0],
             "price": evaluation.price,
             "sequence": evaluation.sequence}
            for evaluation in self.plan
        ], columns=["segment", "offer", "length",
                    "price", "sequence"])
        if sort_by is not None:
            dataframe = dataframe.sort_values(by=sort_by)
        return dataframe

    def to_SeqRecord(self, record=None, record_id="Ordering plan"):
        """
        >>> record = to_SeqRecord(solution)
        >>> # Let's plot with DnaVu:
        >>> from dnavu import create_record_plot
        >>> from bokeh.io import output_file, show
        >>> output_file("view.html")
        >>> plot = create_record_plot(record)
        >>> show(plot)
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("Install Biopython to use to_SeqRecord.")
        offers = sorted(self.plan, key=lambda s: s.segment)
        if record is None:
            record = SeqRecord(Seq(self.full_sequence, DNAAlphabet()),
                               id=record_id)
        else:
            record = deepcopy(record)
        new_features = []
        for offer in offers:
            ind = record.seq.find(offer.sequence)
            if (ind != -1):
                start, end = ind, ind + len(offer.sequence)
            else:
                start, end = offer.segment[0], offer.segment[1]
            feature = SeqFeature(
                FeatureLocation(start, end, 1), type="order",
                qualifiers={
                    "offer": offer.offer, "price": offer.price,
                    "locus_tag": "%s, %d$" % (offer.offer, offer.price)
                }
            )
            new_features.append(feature)
        record.features = new_features + record.features
        return record

    def to_genbank(self, filename, record=None, record_id="Ordering plan"):
        record = self.to_SeqRecord(record=record, record_id=record_id)
        with open(filename, "w+") as f:
            SeqIO.write(record, f, "genbank")

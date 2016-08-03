"""Classes to represent assembly strategies in DnaWeaver"""
import itertools as itt
from copy import deepcopy
from collections import defaultdict

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except:
    PANDAS_AVAILABLE = False

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .optimization import NoSolutionFoundError

class DnaQuote:
    """Class to represent the Quote returned by a DNA source in response to
    a quote request for a given DNA sequence.

    Parameters
    -----------

     source

     sequence

     price

     accepted

     lead_time

     ordering_plan

     metadata

     message

    """

    def __init__(self, source, sequence, price=None, accepted=True,
                 lead_time=None, ordering_plan=None, metadata={}, message=""):
        self.source = source
        self.accepted = accepted
        self.price = price
        self.lead_time = lead_time
        self.sequence = sequence
        self.message = message
        self.ordering_plan = ordering_plan
        self.metadata = metadata

    def __repr__(self):
        if self.accepted:
            infos = "price %.02f" % self.price
            if self.lead_time is not None:
                infos = infos + ", lead_time %0.1f" % self.lead_time
        else:
            infos = "refused: %s" % self.message
        return "From %s, %s" % (self.source, infos)

    def compute_assembly_tree(quote, id_prefix="Op. "):

        counter = itt.count()

        def rec(quote):
            source = quote.source
            if (hasattr(quote.source, "dna_source") or
                    hasattr(quote.source, "primers_dna_source")):
                if quote.ordering_plan is None:
                    quote = source.get_quote(quote.sequence,
                                             max_lead_time=quote.lead_time,
                                             with_ordering_plan=True)
                segments = {
                    segment: rec(subquote)
                    for segment, subquote in
                    sorted(quote.ordering_plan.quotes.items(),
                           key=lambda it: it[0])
                }
                return AssemblyOperation(
                    id="%s%04d" % (id_prefix, counter.next()),
                    quote=quote,
                    segments=segments
                )
            else:
                return AssemblyOperation(
                    id="%s%04d" % (id_prefix, counter.next()),
                    quote=quote,
                    segments={}
                )

        return rec(quote)

class AssemblyOperation:
    """A class to represent assembly operation to model, analyze and report
    assembly trees.

    AssemblyOperations contain "segments" which can originate from other
    AssemblyOperations, so that an AssemblyOperation instance can be the
    root of a tree-like structure.

    Parameters
    ----------

    id
      A string identifying uniquely the assembly operation

    quote
      DnaQuote associated with the operation, indicating price, leadtime,
      sequence, etc.

    segments
      A dict of the form { (start1, end1): op1, (start2, end2): op2, ...}
      where the `(start, end)` represent sub-segments of the final sequence and
      `op1, op2` represent the AssemblyOperations that create the corresponding
      fragments

    deadline
      Deadline for the operation (see the `propagate_deadline` method for this
      class)

    """

    def __init__(self, id, quote, segments, deadline=None):
        self.id = id
        self.quote = quote
        self.segments = segments
        self.deadline = deadline

    def return_tree_as_list(self):
        """Return a list containing the current AssemblyOperation and all its
        sub-operations and their respective sub-operations.

        Said otherwise, it flattens the assembly tree into the list of all
        nodes.
        """
        return [self] + sum([child.return_tree_as_list()
                             for segment, child in self.segments.items()], [])

    def propagate_deadline(self, deadline):
        self.deadline = deadline
        children_deadline = deadline - self.step_duration
        for segment, child in self.segments.items():
            child.propagate_deadline(children_deadline)

    @property
    def step_duration(self):
        if self.segments == {}:
            children_lead_time = 0
        else:
            children_lead_time = max(child.quote.lead_time
                                     for _, child in self.segments.items())
        return self.quote.lead_time - children_lead_time

    def compute_assembly_levels(self):

        levels = defaultdict(lambda: [])
        edges = []

        def rec(subtree, depth=0):

            levels[depth].append(subtree)

            if hasattr(subtree, "segments"):
                for other in subtree.segments.values():
                    edges.append((other, subtree))
                    rec(other, depth + 1)
        rec(self)
        levels = [levels[i] for i in sorted(levels.keys())][::-1]
        levels = [sorted(level, key=lambda e: e.id)[::-1] for level in levels]
        return edges, levels



class DnaOrderingPlan:
    """Class for representing a group of orders and exporting to many formats.

    Parameters
    ----------

    plan
      A dictionary { sequence: offer_evaluation, ... } where `sequence`
      is an ATGC string sequence of DNA, and offer_evaluation is a
      DnaOfferEvaluation.

    """

    def __init__(self, quotes, full_sequence=None):
        self.quotes = quotes
        self.full_sequence = full_sequence

    def summary(self):
        """Return a print-friendly, simple string of the ordering plan."""
        plan = "\n  ".join(str(segment) + ": " + str(quote)
                           for segment, quote in sorted(self.quotes.items()))
        total_price = self.total_price()
        final_txt = "Ordering plan:\n  %s\n  Price:%.02f" % (plan, total_price)
        lead_time = self.overall_lead_time()
        if lead_time is not None:
            final_txt = final_txt + ", lead_time:%.1f" % lead_time
        return final_txt

    def total_price(self):
        return sum(quote.price for quote in self.quotes.values())

    def overall_lead_time(self):
        lead_times = [quote.lead_time for quote in self.quotes.values()]
        if any(lead_time is None for lead_time in lead_times):
            return None
        else:
            return max(lead_times)

    def cuts_indices(self):
        return sorted([offer.segment[1] for offer in self.quotes.keys()])[:-1]

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

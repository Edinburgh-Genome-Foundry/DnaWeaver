from copy import deepcopy
from dnachisel.constraints import Constraint
from dnachisel import DnaCanvas
from .optimization import NoSolutionFoundError
from copy import copy

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


class DnaQuote:

    def __init__(self, source, sequence, price=None, accepted=True, delay=None,
                 ordering_plan=None, message=""):
        self.source = source
        self.accepted = accepted
        self.price = price
        self.delay = delay
        self.sequence = sequence
        self.message = message
        self.ordering_plan = ordering_plan

    def __repr__(self):
        if self.accepted:
            infos = "price %.02f" % self.price
            if self.delay is not None:
                infos = infos + ", delay %0.1f" % self.delay
        else:
            infos = "refused: %s" % self.message
        return "From %s, %s" % (self.source, infos)


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
        delay = self.overall_delay()
        if delay != None:
            final_txt = final_txt + ", Delay:%.1f" % delay
        return final_txt

    def total_price(self):
        return sum(quote.price for quote in self.quotes.values())

    def overall_delay(self):
        delays = [quote.delay for quote in self.quotes.values()]
        if any(delay is None for delay in delays):
            return None
        else:
            return max(delays)

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
from ..DnaQuote import DnaQuote
from .DnaSource import DnaSource

class DnaSourcesComparator(DnaSource):
    """Special source that compares quotes from other DnaSources.

    Upon receiving a sequence, that source will submit the sequence to
    Downstream sources, which deliver each one optimal quote. The comparator
    then returns the one quote with the lowest price.

    Parameters
    ----------

    dna_sources
      List of `DnaSources` that

    memoize
      Whether the quotes should be kept in memory to avoid re-computing the
      same quote several times. Can accelerate computations but is
      RAM-expensive.
    """
    class_description = "DNA sources comparator"
    operation_type = "comparison"
    report_fa_symbol = u""
    report_fa_symbol_plain = u"circle-o"
    report_color = "#000000"

    def __init__(self, dna_sources, memoize=False, sequence_constraints=(),
                 name="comparator"):
        self.dna_sources = dna_sources
        self.memoize = memoize
        self.sequence_constraints = sequence_constraints
        self.memoize_dict = {}
        self.name = name
        self.min_basepair_price = min([source.min_basepair_price
                                       for source in dna_sources])

    def get_best_price(self, sequence, max_lead_time=None,
                       with_assembly_plan=False):
        """Returns a price-optimal DnaQuote for the given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_assembly_plan
          If True, the assembly plan is added to the quote
        """

        best_quote = DnaQuote(self, sequence, accepted=False,
                              message="Sequence was rejected by all sources.")
        best_basepair_price = 1e8
        for source in sorted(self.dna_sources,
                             key=lambda s: s.min_basepair_price):
            # print (self.name, len(sequence), source)
            if source.min_basepair_price > best_basepair_price:
                break
            quote = source.get_quote(sequence, max_lead_time=max_lead_time,
                                     with_assembly_plan=with_assembly_plan)
            # print ('=>', quote.accepted, quote.price)
            if not quote.accepted:
                continue
            quote_basepair_price = quote.price / float(len(sequence))
            if quote_basepair_price < best_basepair_price:
                best_quote = quote
                best_basepair_price = quote_basepair_price
        # print ('---- went for ', best_quote.source)
        return best_quote

    def additional_dict_description(self):
        return {
            "dna sources": [source.name for source in self.dna_sources],
        }

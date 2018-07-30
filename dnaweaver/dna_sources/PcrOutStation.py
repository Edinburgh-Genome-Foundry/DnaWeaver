from ..DnaQuote import DnaQuote
from .DnaSource import DnaSource
from ..constraints import SequenceLengthConstraint
from ..biotools import (reverse_complement, largest_common_substring,
                        blast_sequence)
from ..tools import functions_list_to_string
from .DnaSourcesComparator import DnaSourcesComparator

class PcrOutStation(DnaSource):
    """Class to represent databases of constructs which can be (in part) reused

    A blast database contains the sequences of all available constructs.
    Given a sequence, the PcrOutStation finds whether it is possible to order
    two primers to extract this sequence from the constructs in the BLAST
    database.

    Parameters
    ----------

    name
      Name of the PCR station (e.g. "Lab constructs PCR station")

    primers_dna_source
      DnaSource providing the primers (will typically be an CommercialDnaOffer)

    blast_database

    sequences
      A dictionary {seq_name: sequence_in_atgc}

    pcr_homology_length

    max_overhang_length

    extra_cost

    extra_time

    max_amplicon_length

    blast_word_size

    memoize

    sequence_constraints

    """
    class_description = "PCR-out station"
    operation_type = "PCR"
    report_fa_symbol = u""
    report_fa_symbol_plain = "exchange"
    report_color = "#eeffee"

    def __init__(self, name, primers_dna_source, blast_database=None,
                 sequences=None, pcr_homology_length=25,
                 max_overhang_length=40, extra_cost=0, extra_time=0,
                 max_amplicon_length=None, blast_word_size=50, memoize=False,
                 sequence_constraints=()):
        self.name = name
        self.blast_database = blast_database
        self.set_suppliers(primers_dna_source)
        self.pcr_homology_length = pcr_homology_length
        self.max_overhang_length = max_overhang_length
        self.extra_time = extra_time
        self.extra_cost = extra_cost
        self.max_amplicon_length = max_amplicon_length
        self.blast_word_size = blast_word_size
        self.sequence_constraints = list(sequence_constraints)
        if max_amplicon_length is not None:
            c = SequenceLengthConstraint(max_length=max_amplicon_length)
            self.sequence_constraints = [c] + self.sequence_constraints
        self.memoize = memoize
        self.memoize_dict = {}
        self.sequences = sequences

    def get_hits(self, sequence):
        """Return the hits of the given sequence against the blast database
        (in Biopython format)"""
        if self.sequences is not None:
            result = []
            for dna_name, seq in self.sequences.items():
                match_coords = largest_common_substring(sequence, seq,
                                                        self.max_overhang_length)
                if match_coords:
                    result.append((dna_name, match_coords, None))
            return result
        else:
            record = blast_sequence(sequence, self.blast_database,
                                    perc_identity=99,
                                    use_megablast=True,
                                    word_size=self.blast_word_size)

            return [
                (al.hit_id + "_h%03d" % i,
                 (hit.query_start, hit.query_end),
                 hit.sbjct)
                for al in record.alignments
                for i, hit in enumerate(al.hsps)
            ]

    def get_best_price(self, sequence, max_lead_time=None,
                       with_assembly_plan=False):
        """Return a price-optimal DnaQuote for the given sequence.

        It will find a possible hit in the blast database, find the primers to
        order for the PCR, compute the overall price and lead time, and return
        a quote.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quote

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_assembly_plan
          If True, the assembly plan is added to the quote
        """

        for subject, (hit_start, hit_end), _ in self.get_hits(sequence):

            largest_overhang = max(hit_start, len(sequence) - hit_end)

            if largest_overhang <= self.max_overhang_length:
                primer_l_end = hit_start + self.pcr_homology_length
                primer_left = sequence[:primer_l_end]
                primer_r_end = hit_end - self.pcr_homology_length
                primer_right = reverse_complement(sequence[primer_r_end:])

                primer_max_lead_time = (None if max_lead_time is None else
                                        max_lead_time - self.extra_time)
                quotes = [
                    self.primers_dna_source.get_quote(
                        primer, max_lead_time=primer_max_lead_time
                    )
                    for primer in [primer_left, primer_right]
                ]
                if not all(quote.accepted for quote in quotes):
                    continue  # primers inorderable

                if max_lead_time is not None:
                    overall_lead_time = (max(quote.lead_time
                                             for quote in quotes) +
                                         self.extra_time)
                else:
                    overall_lead_time = None
                total_price = (sum(quote.price for quote in quotes) +
                               self.extra_cost)

                if with_assembly_plan:
                    assembly_plan = {
                        (0, primer_l_end): quotes[0],
                        (primer_r_end, len(sequence)): quotes[1]
                    }
                else:
                    assembly_plan = None

                return DnaQuote(self, sequence, accepted=True,
                                lead_time=overall_lead_time,
                                price=total_price,
                                assembly_plan=assembly_plan,
                                metadata={"subject": subject,
                                          "location": (hit_start, hit_end)})
        return DnaQuote(self, sequence, accepted=False,
                        message="No valid match found")

    def suggest_cuts(self, sequence):
        suggested_cuts = []
        for name, subseq in self.sequences:
            if subseq in sequence:
                index = sequence.find(subseq)
                suggested_cuts += [index, index + len(subseq)]
        return sorted(list(set(suggested_cuts)))

    def pre_blast(self, sequence):
        """Pre-compute the BLAST of the current sequence against the database.

        Once a pre-blast has been performed, this PcrOutStation becomes
        specialized on that sequence and its subsequences, do not feed it with
        another different sequence. Do `self.sequences=None` to reinitialize
        and de-specialize this PcrOutStation.

        Examples
        --------

        >>> pcr_station = PcrOutStation("some_blast_database")
        >>> top_station = # some assembly station depending on pcr_station
        >>> pcr_station.pre-blast(my_sequence)
        >>> top_station.get_quote(my_sequence)
        >>> pcr_station.sequences=None # de-specializes the pcr station.
        """
        self.sequences = None  # destroy current pre-blast (used by get_hits)
        self.sequences = {
            subject: sbjct
            for subject, (start, end), sbjct in self.get_hits(sequence)
        }

    def additional_dict_description(self):
        constraints = functions_list_to_string(self.sequence_constraints)
        return {
            "BLAST database": self.blast_database,
            "primers dna source": self.primers_dna_source.name,
            "pcr homology length": self.pcr_homology_length,
            "max overhang length": self.max_overhang_length,
            "extra time": self.extra_time,
            "extra cost": self.extra_cost,
            "max amplicon length": self.max_amplicon_length,
            "BLAST word size": self.blast_word_size,
            "sequence constraints": constraints,
        }

    @classmethod
    def from_dict(cls, data):
        if 'dna_bank' in data:
            blast_database = cls.dna_banks[data['dna_bank']]
        else:
            blast_database = data['blast_database']
        return PcrOutStation(
            name=data['name'],
            primers_dna_source=data['suppliers'],
            pcr_homology_length=data['pcr_homology_length'],
            max_overhang_length=data['max_overhang_length'],
            max_amplicon_length=data['max_amplicon_length'],
            extra_cost=data['cost'],
            blast_database=blast_database
        )

    def set_suppliers(self, suppliers):
        if hasattr(suppliers, '__iter__'):
            if len(suppliers) > 1:
                self.primers_dna_source = DnaSourcesComparator(
                    name=self.name + ' comparator', suppliers=suppliers)
            else:
                self.primers_dna_source = suppliers[0]
        else:
            self.primers_dna_source = suppliers

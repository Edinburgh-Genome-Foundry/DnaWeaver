from ...DnaQuote import DnaQuote
from ...SegmentSelector import TmSegmentSelector, FixedSizeSegmentSelector
from ...biotools import (
    reverse_complement,
    largest_common_substring,
    blast_sequence,
    perfect_match_locations_in_hsp,
)
from ...tools import functions_list_to_string
from ..DnaSupplier import DnaSupplier
from ..builtin_constraints import SequenceLengthConstraint
from .DnaSuppliersComparator import DnaSuppliersComparator


class PcrExtractionStation(DnaSupplier):
    """Class to represent databases of constructs which can be (in part) reused

    A blast database contains the sequences of all available constructs.
    Given a sequence, the PcrExtractionStation finds whether it is possible to order
    two primers to extract this sequence from the constructs in the BLAST
    database.

    Parameters
    ----------

    name
      Name of the PCR station (e.g. "Lab constructs PCR station").

    primers_supplier
      DnaSupplier providing the primers (will typically be an CommercialDnaOffer).

    blast_database

    sequences
      A dictionary {seq_name: sequence_in_atgc}.

    pcr_homology_length

    max_overhang_length

    extra_cost

    extra_time

    max_amplicon_length

    blast_word_size
      BLAST parameter, the larger number will lead to faster BLASTing but with
      a risk of missing short sequences of length shorter than this value.

    memoize

    sequence_constraints

    """

    class_description = "PCR-extraction station"
    operation_type = "PCR"
    report_fa_symbol = u""
    report_fa_symbol_plain = "recycle"
    report_color = "#eeffee"
    dna_banks = {}

    def __init__(
        self,
        name,
        primers_supplier,
        homology_selector,
        blast_database=None,
        sequences=None,
        pcr_homology_length=25,
        max_overhang_length=40,
        extra_cost=0,
        extra_time=0,
        max_amplicon_length=None,
        blast_word_size=50,
        memoize=False,
        sequence_constraints=(),
    ):
        self.name = name
        self.homology_selector = homology_selector
        self.blast_database = blast_database
        self.set_suppliers(primers_supplier)
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
            self.min_basepair_price = (
                2 * 20 * self.primers_supplier.min_basepair_price
            ) / self.max_amplicon_length
        self.memoize = memoize
        self.memoize_dict = {}
        self.sequences = sequences

    def _get_hits(self, sequence):
        """Return the hits of the given sequence against the blast database
        in format [(part_hit, (start, end), subsequence).
        """
        if self.sequences is not None:
            result = []
            for dna_name, seq in self.sequences.items():
                match_coords = largest_common_substring(
                    sequence, seq, self.max_overhang_length
                )
                # if match_coords
                if match_coords:
                    result.append((dna_name, match_coords, None))
            return result
        else:
            record = blast_sequence(
                sequence,
                self.blast_database,
                perc_identity=98,
                use_megablast=True,
                word_size=self.blast_word_size,
            )

            return [
                (al.hit_id + "_h%03d" % i, (hit.query_start, hit.query_end), hit.sbjct,)
                for al in record.alignments
                for i, hit in enumerate(al.hsps)
            ]

    def blast_sequence(self, sequence, cutoff=40):
        if hasattr(sequence, "seq"):
            sequence = str(sequence.seq)
        record = blast_sequence(
            sequence,
            self.blast_database,
            perc_identity=98,
            use_megablast=True,
            word_size=self.blast_word_size,
        )
        return [
            (al.hit_id + "_h%03d_%02d" % (i, j), (start, end), sequence[start:end],)
            for al in record.alignments
            for i, hit in enumerate(al.hsps)
            for j, (start, end) in enumerate(
                perfect_match_locations_in_hsp(hit, span_cutoff=cutoff)
            )
        ]

    def get_best_price(
        self, sequence, max_lead_time=None, with_assembly_plan=False,
    ):
        """Return a price-optimal DnaQuote for the given sequence.

        It will find a possible hit in the blast database, find the primers to
        order for the PCR, compute the overall price and lead time, and return
        a quote.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quote.

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        with_assembly_plan
          If True, the assembly plan is added to the quote.
        """
        hits = self._get_hits(sequence)

        # For each hit in the database, see if there is a
        for subject, (hit_start, hit_end), _ in hits:
            if min(hit_start, hit_end) < 0:
                continue
            largest_overhang = max(hit_start, len(sequence) - hit_end)

            if largest_overhang > self.max_overhang_length:
                continue
            for i in range(min(len(sequence), self.max_overhang_length) - hit_start):
                subseq = sequence[hit_start + i :]
                # print (hit_start, i, subseq, hit_end, hit_start)
                left_location = self.homology_selector.compute_segment_location(
                    subseq, 0
                )
                if left_location is not None:
                    l_start, l_end = left_location
                    primer_left = sequence[: l_end + hit_start + i]
                    break
            else:
                continue
            right_padding = len(sequence) - hit_end
            for i in range(self.max_overhang_length - right_padding):
                subseq = sequence[: hit_end - i]
                right_location = self.homology_selector.compute_segment_location(
                    subseq, hit_end - i
                )
                if right_location is not None:
                    r_start, r_end = right_location
                    primer_right = reverse_complement(sequence[r_start:])
                    break
            else:
                continue
            # primer_l_end = hit_start + self.pcr_homology_length
            # primer_left = sequence[:primer_l_end]
            # primer_r_end = hit_end - self.pcr_homology_length
            # primer_right = reverse_complement(sequence[primer_r_end:])

            primer_max_lead_time = (
                None if max_lead_time is None else max_lead_time - self.extra_time
            )
            quotes = [
                self.primers_supplier.get_quote(
                    primer, max_lead_time=primer_max_lead_time
                )
                for primer in [primer_left, primer_right]
            ]
            if not all(quote.accepted for quote in quotes):
                continue  # primers inorderable

            if max_lead_time is not None:
                overall_lead_time = (
                    max(quote.lead_time for quote in quotes) + self.extra_time
                )
            else:
                overall_lead_time = None
            total_price = sum(quote.price for quote in quotes) + self.extra_cost

            if with_assembly_plan:
                assembly_plan = {
                    (0, int(l_end)): quotes[0],
                    (int(r_start), len(sequence)): quotes[1],
                }
            else:
                assembly_plan = None

            return DnaQuote(
                self,
                sequence,
                accepted=True,
                lead_time=overall_lead_time,
                price=total_price,
                assembly_plan=assembly_plan,
                message="From %s" % subject,
                metadata={"subject": subject, "location": (hit_start, hit_end),},
            )
        if len(hits):
            message = "Some matches found but could not find suitable primer design."
        else:
            message = "No BLAST hit found"
        return DnaQuote(self, sequence, accepted=False, message=message)

    def suggest_cuts(self, sequence):
        if self.sequences is None:
            return []
        suggested_cuts = []
        for name, subseq in self.sequences.items():
            index = sequence.find(subseq)
            if index >= 0:
                suggested_cuts.extend([index, index + len(subseq)])
        return sorted(list(set(suggested_cuts)))

    def prepare_on_sequence(self, sequence):
        """Pre-compute the BLAST of the current sequence against the database.

        Once a pre-blast has been performed, this PcrExtractionStation becomes
        specialized on that sequence and its subsequences, do not feed it with
        another different sequence. Do `self.sequences=None` to reinitialize
        and de-specialize this PcrExtractionStation.

        Examples
        --------

        >>> pcr_station = PcrExtractionStation("some_blast_database")
        >>> top_station = # some assembly station depending on pcr_station
        >>> pcr_station.prepare_on_sequence(my_sequence)
        >>> top_station.get_quote(my_sequence)
        >>> pcr_station.sequences=None # de-specializes the pcr station.
        """
        if self.blast_database is None:
            return
        hits = self.blast_sequence(sequence)
        self.sequences = None  # destroy current pre-blast (used by _get_hits)
        self.sequences = {subject: seq for subject, (start, end), seq in hits}

    def additional_dict_description(self):
        constraints = functions_list_to_string(self.sequence_constraints)
        return {
            "BLAST database": self.blast_database,
            "primers dna source": self.primers_supplier.name,
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
        if "dna_bank" in data:
            blast_database = cls.dna_banks[data["dna_bank"]]
        else:
            blast_database = data["blast_database"]

        if data["homology_type"] == "tm":
            min_oh_size, max_oh_size = data["homology_size_range"]
            min_tm, max_tm = data["tm_range"]
            homology_selector = TmSegmentSelector(
                min_size=min_oh_size,
                max_size=max_oh_size,
                min_tm=min_tm,
                max_tm=max_tm,
            )
        else:
            homology_selector = FixedSizeSegmentSelector(
                segment_size=data["homology_size"]
            )
        return PcrExtractionStation(
            name=data["name"],
            homology_selector=homology_selector,
            primers_supplier=data["suppliers"],
            max_overhang_length=data["max_overhang_length"],
            max_amplicon_length=data["max_amplicon_length"],
            extra_cost=data["cost"],
            blast_database=blast_database,
        )

    def set_suppliers(self, suppliers):
        if hasattr(suppliers, "__iter__"):
            if len(suppliers) > 1:
                self.primers_supplier = DnaSuppliersComparator(
                    name=self.name + " comparator", suppliers=suppliers
                )
            else:
                self.primers_supplier = suppliers[0]
        else:
            self.primers_supplier = suppliers
        self.suppliers = [self.primers_supplier]  # for network reconstitution

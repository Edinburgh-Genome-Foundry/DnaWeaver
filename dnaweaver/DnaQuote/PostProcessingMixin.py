import itertools as itt
import os
import tempfile
from ..biotools import blast_sequence


class PostProcessingMixin:
    def compute_full_assembly_plan(self, id_prefix="S", id_digits=5):
        """ """

        counter = itt.count()

        def rec(quote):
            if not quote.accepted:
                return quote

            if any(
                [
                    hasattr(quote.source, attr)
                    for attr in ["supplier", "primers_supplier"]
                ]
            ):
                if quote.assembly_plan is None:
                    quote = quote.source.get_quote(
                        quote.sequence,
                        max_lead_time=quote.lead_time,
                        with_assembly_plan=True,
                    )
                segments = {
                    segment: rec(subquote)
                    for segment, subquote in sorted(
                        quote.assembly_plan.items(), key=lambda item: item[0]
                    )
                }
                quote.assembly_plan = segments
            if id_prefix:
                index = next(counter)
                quote.id = "{id_prefix}_{index:0{id_digits}}".format(
                    id_prefix=id_prefix, index=index, id_digits=id_digits
                )
            return quote

        rec(self)

        if id_prefix:
            index = next(counter)
            self.id = "{id_prefix}_{index:0{id_digits}}".format(
                id_prefix=id_prefix, index=index, id_digits=id_digits
            )
        self.full_assembly_plan_computed = True

    def compute_fragments_final_locations(self):
        """Compute the exact final location of the fragments in the final
        sequence.
        """
        if not self.full_assembly_plan_computed:
            self.compute_full_assembly_plan()
        quotes = self.tree_as_list()
        quotes_dict = {quote.id: quote for quote in quotes}
        _, temp_fasta = tempfile.mkstemp(suffix=".fa")
        with open(temp_fasta, "w+") as f:
            for quote in quotes:
                f.write(">%s\n%s\n" % (quote.id, quote.sequence))
        results = blast_sequence(
            self.sequence, subject=temp_fasta, word_size=10, perc_identity=100
        )

        if isinstance(results, list):
            alignments = sum([rec.alignments for rec in results], [])
        else:
            alignments = results.alignments

        for al in alignments:
            hit = max(al.hsps, key=lambda hit: hit.align_length)
            final_location = sorted((hit.query_start, hit.query_end))
            matching_segment = sorted((hit.sbjct_start, hit.sbjct_end))
            quotes_dict[al.hit_def].final_location = final_location
            quotes_dict[al.hit_def].matching_segment = matching_segment
        os.remove(temp_fasta)

    def propagate_deadline(self, deadline):
        """Add a `deadline` attribute to the quote and propagate it to
        the quote's children by taking into account the duration of operations.

        For instance if "self" has a duration of 5 and receives a deadline
        of 8, the quotes that "self" depends on will receive a deadline of
        8-5=3.
        """
        self.deadline = deadline
        children_deadline = deadline - self.step_duration
        if self.assembly_plan is not None:
            for segment, child in self.assembly_plan.items():
                child.propagate_deadline(children_deadline)

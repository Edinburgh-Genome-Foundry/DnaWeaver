"""Classes to represent assembly strategies in DnaWeaver."""


from collections import defaultdict
from .ExportsMixin import ExportsMixin
from .PostProcessingMixin import PostProcessingMixin


class DnaQuote(ExportsMixin, PostProcessingMixin):
    """Class to represent the Quote returned by a DNA source in response to
    a quote request for a given DNA sequence.

    Parameters
    ----------

    source
      The DnaSupplier which issued the quote.

    sequence
      The sequence being quoted.

    accepted
      False if the sequence was rejected by the supplier, else True.

    price
      Amount asked by the DnaSupplier to build the requested sequence.

    lead_time
      Number of days required by the DnaSupplier to deliver the sequence.

    assembly_plan
      A dict of the form { (start1, end1): op1, (start2, end2): op2, ...}
      where the `(start, end)` represent sub-segments of the final sequence and
      `op1, op2` represent the AssemblyOperations that create the corresponding
      fragments.

    message
      A message added to the quote for more information, for instance the
      reason of the rejection.

    id
      A string identifying this assembly operation.

    deadline
      If the assembly plan was computed with a dealine in mind, this is the
      number of days at which the sequence should be delivered to ensure that
      the project will be finished on time (so the global deadline minus the
      duration of downstream operations). This is for drawing assembly
      timelines. (see the `propagate_deadline` method for this class).

    metadata
      Any other data on the quote.
    """

    def __init__(
        self,
        source,
        sequence,
        price=None,
        accepted=True,
        lead_time=None,
        assembly_plan=None,
        metadata=None,
        message="",
        id=None,
        deadline=None,
    ):
        self.source = source
        self.accepted = accepted
        self.price = price
        self.lead_time = lead_time
        self.sequence = sequence
        self.message = message
        self.assembly_plan = assembly_plan
        self.full_assembly_plan_computed = False
        self.metadata = metadata or {}
        self.id = id
        self.deadline = None

    def __repr__(self):
        if self.accepted:
            infos = "price %.02f" % self.price
            if self.lead_time is not None:
                infos = infos + " - lead_time %0.1f" % self.lead_time
        else:
            infos = "refused: %s" % self.message
        result = "From %s - %s" % (self.source, infos)
        if len(self.message):
            result += " - " + self.message
        return result

    @property
    def step_duration(self):
        """Duration of the final step of the assembly in the quote.

        Inferred by substracting the lead time of children operations to the
        lead time of this operation.
        """
        children_lead_time = self.children_overall_lead_time()
        if children_lead_time is None:
            return self.lead_time
        else:
            return self.lead_time - self.children_overall_lead_time()

    def compute_assembly_levels(self):
        """Return edges and levels for drawing the assembly graph nicely."""

        levels = defaultdict(lambda: [])
        edges = []

        def rec(subtree, depth=0):

            levels[depth].append(subtree)

            if subtree.assembly_plan is not None:
                for other in subtree.assembly_plan.values():
                    edges.append((other, subtree))
                    rec(other, depth + 1)

        rec(self)
        levels = [levels[i] for i in sorted(levels.keys())][::-1]
        levels = [sorted(level, key=lambda e: e.id)[::-1] for level in levels]
        return edges, levels

    def assembly_step_summary(self):
        """Return a print-friendly, simple string of the ordering plan."""
        plan = "\n  ".join(
            "%s-%s: %s" % (start, end, str(quote))
            for (start, end), quote in sorted(self.assembly_plan.items())
        )
        title = "Ordering plan (%s):" % self.source.name
        final_txt = "%s:\n  %s\nPrice:%.02f" % (title, plan, self.price)
        if self.lead_time is not None:
            final_txt = final_txt + ", total lead_time:%.1f" % self.lead_time
        return final_txt

    def children_total_price(self):
        """Return the total price of all sub-operations (apart from current)."""
        # print ([quote.accepted for quote in self.assembly_plan.values()])
        return sum(quote.price for quote in self.assembly_plan.values())

    def children_overall_lead_time(self):
        """Return the max lead time of all sub-operation (current one not
        included)."""
        if self.assembly_plan is None:
            return None
        lead_times = [quote.lead_time for quote in self.assembly_plan.values()]
        if len(lead_times) == 0:
            return None
        elif any(lead_time is None for lead_time in lead_times):
            return None
        else:
            return max(lead_times)

    def cuts_indices(self):
        """Return the locations of the `cuts` where the sequence is assembled.
        """
        return sorted([segment[1] for segment in self.assembly_plan])[:-1]

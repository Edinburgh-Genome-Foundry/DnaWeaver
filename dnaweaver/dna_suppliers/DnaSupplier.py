# -*- coding: utf-8 -*-

from collections import defaultdict

import numpy as np
from ..DnaQuote import DnaQuote

# Attempt import of optional module DNA Chisel
try:
    from dnachisel import Specification, DnaOptimizationProblem

    DNACHISEL_AVAILABLE = True
except:
    Specification = type(
        "BLANK"
    )  # meant to be a fake type that matches nothing
    DNACHISEL_AVAILABLE = False


class DnaSupplier:
    """Base class for all DnaSuppliers, which are the elements of the supply
    networks used to define assembly problems in DnaWeaver."""

    min_basepair_price = 0

    def get_quote(
        self,
        sequence,
        max_lead_time=None,
        max_price=None,
        with_assembly_plan=False,
        time_resolution=1.0,
    ):
        """Return a DnaQuote with price, lead time, etc. for a given sequence.

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times
          Note that either this parameter or `max_lead_time` must be None

        with_assembly_plan
          If True, the assembly plan is added to the quote

        time_resolution
          Time resolution for the bisecting search if `max_price` is not None.

        Returns
        -------

        A DnaQuote object.

        """

        if max_lead_time is None:
            max_lead_time = np.inf

        if self.memoize:
            args = (
                hash(sequence),
                max_lead_time,
                max_price,
                with_assembly_plan,
            )
            quote = self.memoize_dict.get(args, None)
            if quote is not None:
                return quote

        if not self.verify_constraints(sequence):
            quote = DnaQuote(
                self,
                sequence,
                accepted=False,
                message="Sequence does not pass constraints.",
            )

        elif max_price is not None:
            quote = self.get_best_lead_time_under_price_limit(
                sequence,
                max_price=max_price,
                time_resolution=time_resolution,
                with_assembly_plan=with_assembly_plan,
            )
        else:
            quote = self.get_best_price(
                sequence,
                max_lead_time=max_lead_time,
                with_assembly_plan=with_assembly_plan,
            )

        if self.memoize:
            self.memoize_dict[args] = quote

        return quote

    def __repr__(self):
        return self.name

    def get_best_lead_time_under_price_limit(
        self, sequence, max_price, time_resolution, with_assembly_plan=False
    ):
        """Return the quote with fastest lead time under the budget constraint

        Parameters
        ----------



        sequence (str)
          The sequence submitted to the Dna Source for a quots

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times
          Note that either this parameter or `max_lead_time` must be None

        with_assembly_plan
          If True, the assembly plan is added to the quote

        time_resolution
          Time resolution for the bisecting search if `max_price` is not None.

        """

        def f(max_lead_time):
            return self.get_quote(
                sequence,
                max_lead_time=max_lead_time,
                with_assembly_plan=with_assembly_plan,
            )

        quote = f(None)
        if (not quote.accepted) or (quote.price > max_price):
            return DnaQuote(
                self,
                sequence,
                accepted=False,
                message="Price over limit even without time limit.",
            )
        step_size = quote.lead_time / 2.0
        lead_time = quote.lead_time - step_size
        best_quote = quote
        while step_size > time_resolution:
            quote = f(lead_time)
            if quote.accepted and (quote.price <= max_price):
                best_quote = quote
                lead_time -= step_size
            else:
                lead_time += step_size
            step_size /= 2.0

        return best_quote

    def verify_constraints(self, sequence):
        """Return True iff `sequence` passes all `self.sequence_constraints`

        Will automatically process DNA-Chisel constraints that would be in
        `self.sequence_constraints`

        """
        constraints = self.sequence_constraints
        if not hasattr(self, "dnachisel_constraints"):
            self.dnachisel_constraints = [
                constraint
                for constraint in self.sequence_constraints
                if isinstance(constraint, Specification)
            ]

        if self.dnachisel_constraints != []:
            if not DNACHISEL_AVAILABLE:
                raise ImportError(
                    "Spotted DNA Chisel constraints, while "
                    "DNA Chisel is not installed."
                )
            # We provide an empty mutation space so it won't be recomputed
            # (which would take time and is useless here!)
            problem = DnaOptimizationProblem(
                sequence, self.dnachisel_constraints, mutation_space=[]
            )
            constraints = [
                constraint
                for constraint in constraints
                if not isinstance(constraint, Specification)
            ] + [lambda seq: problem.all_constraints_pass()]

        return all(constraint(sequence) for constraint in constraints)

    def compute_supply_graph(self):
        """Return elements to plot the supply graph underlying this DnaSupplier

        Returns
        -------

        edges
          A list [(s1,s2), (s1,s3), (s2, s5)...] of couples of DnaSuppliers in a
          supplier-supplied relationship.

        levels
          A list of lists [[s1,s2], [s4,s8,s9]...] of sources. The first
          sublist (first level) are all sources at the farthest distance from
          the current source in the supply graph, and the last sublist contains
          only the current DnaSupplier.
        """

        source_max_level = {}
        edges = []

        def rec(source, depth, seen_sources):
            if source in seen_sources:
                return
            if source not in source_max_level:
                source_max_level[source] = depth
            else:
                source_max_level[source] = max(source_max_level[source], depth)
            new_seen_sources = seen_sources + [source]
            if hasattr(source, "suppliers"):
                for other in source.suppliers:
                    edges.append((other, source))
                    rec(other, depth + 1, new_seen_sources)
            elif hasattr(source, "supplier"):
                edges.append((source.supplier, source))
                rec(source.supplier, depth + 1, new_seen_sources)
            if hasattr(source, "primers_supplier"):
                edges.append((source.primers_supplier, source))
                rec(source.primers_supplier, depth + 1, new_seen_sources)
        
        rec(self, depth=0, seen_sources=[])
        levels = [
            [
                source
                for source, level in source_max_level.items()
                if level == i
            ]
            for i in range(max(source_max_level.values()) + 1)
        ][::-1]

        # seen_sources = set()
        # levels = defaultdict(lambda: [])
        # edges = []

        # def rec(source, depth=0):

        #     if source in seen_sources:
        #         return
        #     seen_sources.add(source)

        #     levels[depth].append(source)
        #     if hasattr(source, "suppliers"):
        #         for other in source.suppliers:
        #             edges.append((other, source))
        #             rec(other, depth + 1)
        #     elif hasattr(source, "supplier"):
        #         edges.append((source.supplier, source))
        #         rec(source.supplier, depth + 1)
        #     if hasattr(source, "primers_supplier"):
        #         edges.append((source.primers_supplier, source))
        #         rec(source.primers_supplier, depth + 1)

        # rec(self)

        # levels = [levels[i] for i in sorted(levels.keys())][::-1]
        return edges, levels

        # seen_sources = set()
        # levels = defaultdict(lambda: [])
        # edges = []

        # def rec(source, depth=0):

        #     if source in seen_sources:
        #         return
        #     seen_sources.add(source)

        #     levels[depth].append(source)
        #     if hasattr(source, "suppliers"):
        #         for other in source.suppliers:
        #             edges.append((other, source))
        #             rec(other, depth + 1)
        #     elif hasattr(source, "supplier"):
        #         edges.append((source.supplier, source))
        #         rec(source.supplier, depth + 1)
        #     if hasattr(source, "primers_supplier"):
        #         edges.append((source.primers_supplier, source))
        #         rec(source.primers_supplier, depth + 1)

        # rec(self)

        # levels = [levels[i] for i in sorted(levels.keys())][::-1]
        # return edges, levels

    def prepare_network_on_sequence(self, sequence):
        edges, levels = self.compute_supply_graph()
        for level in levels:
            for source in level:
                if hasattr(source, "prepare_on_sequence"):
                    source.prepare_on_sequence(sequence)

    def dict_supply_graph(self):
        sources = {}

        def rec(source, depth=0):

            if source in sources:
                return
            if hasattr(source, "is_ghost_source"):
                return
            sources[source.name] = source.dict_description()
            sources[source.name]["_depth"] = depth

            providers = sources[source.name]["providers"] = []
            if hasattr(source, "suppliers"):
                for other in source.suppliers:
                    providers.append(other.name)
                    rec(other, depth + 1)
            if hasattr(source, "dna_supplier"):
                providers.append(source.dna_supplier.name)
                rec(source.dna_supplier, depth + 1)
            if hasattr(source, "primers_supplier"):
                providers.append(source.primers_supplier.name)
                rec(source.primers_supplier, depth + 1)
            if hasattr(source, "dna_suppliers"):
                for other in source.dna_suppliers:
                    providers.append(other.name)
                    rec(other, depth + 1)

        rec(self)
        return sources

    def dict_description(self):
        result = {
            "name": self.name,
            "operation_type": self.operation_type,
            "class": self.class_description,
            "_report_color": self.report_color,
            "_report_fa_symbol": self.report_fa_symbol,
            "_report_fa_symbol_plain": self.report_fa_symbol_plain,
        }
        result.update(self.additional_dict_description())
        return result

    def additional_dict_description(self):
        return {}

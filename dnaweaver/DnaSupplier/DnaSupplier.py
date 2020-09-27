# -*- coding: utf-8 -*-

from collections import defaultdict

import numpy as np
from ..DnaQuote import DnaQuote
from ..biotools import SequenceString

from . import mixins


class DnaSupplier(
    mixins.JsonImportMixin, mixins.ConstraintsMixin, mixins.SupplyNetworkMixin
):
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
          The sequence submitted to the Dna Source for a quote.

        max_lead_time (float)
          If provided, the quote returned is the best quote (price-wise) whose
          lead time is less or equal to max_lead_time.

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times.
          Note that either this parameter or `max_lead_time` must be None.

        with_assembly_plan
          If True, the assembly plan is added to the quote.

        time_resolution
          Time resolution for the bisecting search if `max_price` is not None.

        Returns
        -------

        A DnaQuote object.
        """
        if hasattr(sequence, "seq"):
            # Use a simpler format which still remembers the topology
            sequence = SequenceString.from_record(sequence)

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
        self, sequence, max_price, time_resolution, with_assembly_plan=False,
    ):
        """Return the quote with fastest lead time under the budget constraint/

        Parameters
        ----------

        sequence (str)
          The sequence submitted to the Dna Source for a quote.

        max_price (float)
          If provided, the quote returned is the least-lead-time quote
          whose price is below or equal to `max_price`.
          This is done using bisection and can be slow as it requires to
          re-compute the problem many times.
          Note that either this parameter or `max_lead_time` must be None.

        with_assembly_plan
          If True, the assembly plan is added to the quote.

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

    def prepare_network_on_sequence(self, sequence):
        _edges, levels = self.compute_supply_graph()
        for level in levels:
            for source in level:
                if hasattr(source, "prepare_on_sequence"):
                    source.prepare_on_sequence(sequence)

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

import numpy as np
from ....DnaQuote import DnaQuote
from ....SegmentSelector import TmSegmentSelector, FixedSizeSegmentSelector
from ....DnaAssemblyMethod import (
    OligoAssemblyMethod,
    GibsonAssemblyMethod,
    GoldenGateAssemblyMethod,
)
from ...builtin_constraints import SequenceLengthConstraint
from ...DnaSupplier import DnaSupplier
from ..DnaSuppliersComparator import DnaSuppliersComparator
from .SequenceDecomposer import SequenceDecomposer, NoSolutionFoundError


class DnaAssemblyStation(DnaSupplier):
    """DNA Assembly stations assemble together DNA fragments using a specific
    assembly method.

    Parameters
    ----------

    name
      Name of the station (appears on reports).

    assembly_method
      A DnaAssemblyMethod object specifying how the fragments are assembled,
      what sequences can be assembled, what fragments can be used, etc.

    supplier

    """

    class_description = "DNA assembly station"
    operation_type = "assembly"
    report_fa_symbol = u""
    report_fa_symbol_plain = "flask"
    report_color = "#eeeeff"

    def __init__(
        self,
        name,
        assembly_method,
        supplier,
        memoize=False,
        decomposer_class=None,
        a_star_auto_multiplier=2,
        **solver_kwargs
    ):
        self.name = name
        self.assembly_method = assembly_method
        self.set_suppliers(supplier)
        self.extra_time = assembly_method.duration
        self.extra_cost = assembly_method.cost
        self.sequence_constraints = assembly_method.sequence_constraints
        self.cuts_set_constraints = assembly_method.cuts_set_constraints
        self.memoize = memoize
        if decomposer_class is None:
            self.decomposer_class = SequenceDecomposer
        else:
            self.decomposer_class = decomposer_class
        self.memoize_dict = {}
        self.min_basepair_price = self.supplier.min_basepair_price
        if solver_kwargs.get("a_star_factor", None) == "auto":
            solver_kwargs["a_star_factor"] = (
                a_star_auto_multiplier * self.min_basepair_price
            )
        self.solver_kwargs = solver_kwargs

    def get_quote_for_sequence_segment(
        self, sequence, segment, max_lead_time=None, **kwargs
    ):
        """Return the cost of the segment.

        Is used as the "cost" function for a segment during decomposition
        optimization.

        """
        fragment_to_order = self.assembly_method.compute_fragment_for_sequence_segment(
            sequence=sequence, segment=segment, **kwargs
        )
        return self.supplier.get_quote(fragment_to_order, max_lead_time=max_lead_time)

    def get_assembly_plan_from_cuts(self, sequence, cuts, max_lead_time=None):
        """Return a plan {segment: quote, ...} based on the cut positions.

        Where "segment" is of the form (start, end).
        """
        cuts = sorted(cuts)
        return {
            segment: self.get_quote_for_sequence_segment(
                sequence, segment, max_lead_time=max_lead_time, segment_position=i,
            )
            for i, segment in enumerate(zip(cuts, cuts[1:]))
        }

    def new_sequence_decomposer(
        self,
        sequence,
        max_lead_time=None,
        coarse_grain=1,
        fine_grain=1,
        cut_spread_radius=0,
        a_star_factor=0,
        logger=None,
    ):
        def segment_score(segment):
            quote = self.get_quote_for_sequence_segment(
                sequence, segment, max_lead_time=max_lead_time,
            )
            if not quote.accepted:
                return -1
            else:
                return quote.price

        assembly = self.assembly_method

        return self.decomposer_class(
            sequence_length=len(sequence),
            segment_score_function=segment_score,
            cut_location_constraints=[
                cs(sequence) for cs in assembly.cut_location_constraints
            ],
            segment_constraints=[cs(sequence) for cs in assembly.segment_constraints],
            min_segment_length=assembly.min_segment_length,
            max_segment_length=assembly.max_segment_length,
            forced_cuts=assembly.force_cuts(sequence),
            suggested_cuts=self.compute_suggested_cuts(sequence),
            coarse_grain=coarse_grain,
            fine_grain=fine_grain,
            logger=logger,
            bar_prefix="[%s] " % self.name,
            a_star_factor=a_star_factor,
            path_size_limit=assembly.max_fragments,
            cuts_set_constraints=[cs(sequence) for cs in assembly.cuts_set_constraints],
            cut_spread_radius=cut_spread_radius,
        )

    def compute_suggested_cuts(self, sequence):
        all_supplier_cuts = []
        for supplier in self.suppliers:
            if hasattr(supplier, "suggest_cuts"):
                supplier_cuts = supplier.suggest_cuts(sequence)
                all_supplier_cuts.extend(list(supplier_cuts))
        assembly_method_cuts = self.assembly_method.suggest_cuts(sequence)
        suggested = list(assembly_method_cuts) + list(all_supplier_cuts)
        return sorted(set(suggested))

    def suggest_cuts(self, sequence):
        return self.compute_suggested_cuts(sequence=sequence)

    def get_assembly_plan_for_sequence(
        self,
        sequence,
        max_lead_time=None,
        coarse_grain=1,
        fine_grain=1,
        cut_spread_radius=0,
        a_star_factor=0,
        logger=None,
    ):
        """Return the plan {(seg, ment): quote, ...} of the optimal strategy
        for the sequence's decomposition."""
        # extended_sequence = self.assembly_method.extend_sequence(sequence)
        decomposer = self.new_sequence_decomposer(
            sequence=sequence,
            max_lead_time=max_lead_time,
            coarse_grain=coarse_grain,
            fine_grain=fine_grain,
            a_star_factor=a_star_factor,
            cut_spread_radius=cut_spread_radius,
            logger=logger,
        )
        self._lattest_decomposer = decomposer  # for debugging
        best_cuts = decomposer.compute_optimal_cuts()
        return self.get_assembly_plan_from_cuts(
            sequence, best_cuts, max_lead_time=max_lead_time,
        )

    def get_best_price(
        self, sequence, max_lead_time=None, with_assembly_plan=False,
    ):
        """Returns a price-optimal DnaQuote for the given sequence.

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

        if max_lead_time is not None:
            max_lead_time = max_lead_time - self.extra_time
            if max_lead_time < 0:
                return DnaQuote(
                    self, sequence, accepted=False, message="Lead time limit too short",
                )
        try:
            assembly_plan = self.get_assembly_plan_for_sequence(
                sequence, max_lead_time=max_lead_time, **self.solver_kwargs
            )
        except NoSolutionFoundError as err:
            message = "No solution found! %s" % err
            return DnaQuote(self, sequence, accepted=False, message=message)

        # A solution has been found ! Now compute overall time and lead time.

        if any([q.price is None for q in assembly_plan.values()]):
            return DnaQuote(
                self, sequence, accepted=False, message="No solution found !"
            )
        quote = DnaQuote(self, sequence, assembly_plan=assembly_plan)
        quote.price = quote.children_total_price() + self.extra_cost
        children_lead_time = quote.children_overall_lead_time()
        if (children_lead_time is None) or (self.extra_time is None):
            quote.lead_time = None
        else:
            quote.lead_time = children_lead_time + self.extra_time
        if not with_assembly_plan:
            self.assembly_plan = None
        return quote

    def additional_dict_description(self):
        result = {
            "dna source": self.supplier.name,
            "solver parameters": self.solver_kwargs,
        }
        result.update(
            {
                ("assembly method %s" % k): v
                for k, v in self.assembly_method.dict_description().items()
            }
        )
        return result

    @staticmethod
    def from_dict(data):
        sequence_constraints = []
        if data["use_size_range"]:
            mini, maxi = data["size_range"]
            sequence_constraints.append(SequenceLengthConstraint(mini, maxi))
        min_length, max_length = data["fragments_size_range"]
        if data["method"] == "type_iis":
            gc_range = data.get("overhang_gc_range", [0, 1])
            method = GoldenGateAssemblyMethod(
                duration=data["duration"],
                cost=data["cost"],
                min_overhangs_gc=gc_range[0],
                max_overhangs_gc=gc_range[1],
                enzyme=data["enzyme"],
                max_segment_length=max_length,
                min_segment_length=min_length,
                max_fragments=data["max_fragments"],
                sequence_constraints=sequence_constraints,
            )
        elif data["method"] in [
            "gibson_assembly",
            "yeast_recombination",
            "oligo_assembly",
        ]:
            if data["overhang_type"] == "tm":
                min_oh_size, max_oh_size = data["overhang_size_range"]
                min_tm, max_tm = data["tm_range"]
                overhang_selector = TmSegmentSelector(
                    min_size=min_oh_size,
                    max_size=max_oh_size,
                    min_tm=min_tm,
                    max_tm=max_tm,
                )
            else:
                overhang_selector = FixedSizeSegmentSelector(
                    segment_size=data["overlap"]
                )
            if data["method"] == "oligo_assembly":
                method_class = OligoAssemblyMethod
            else:
                method_class = GibsonAssemblyMethod
            method = method_class(
                duration=data["duration"],
                cost=data["cost"],
                overhang_selector=overhang_selector,
                max_segment_length=max_length,
                min_segment_length=min_length,
                max_fragments=data["max_fragments"],
                sequence_constraints=sequence_constraints,
            )
        if data["use_astar"]:
            a_star_factor = data.get("astar_factor", "auto")
        else:
            a_star_factor = 0
        max_construct_length = max_length * data["max_fragments"]
        if data["grain_type"] == "auto":
            data["coarse_grain"] = int(max_construct_length / 15)
            data["fine_grain"] = int(np.sqrt(data["coarse_grain"]))

        return DnaAssemblyStation(
            name=data["name"],
            supplier=data["suppliers"],
            assembly_method=method,
            coarse_grain=data["coarse_grain"],
            fine_grain=data["fine_grain"],
            a_star_factor=a_star_factor,
            memoize=data.get("memoize", False),
            logger=data.get("logger", None),
        )

    def set_suppliers(self, suppliers):
        if hasattr(suppliers, "__iter__"):
            if len(suppliers) > 1:
                self.supplier = DnaSuppliersComparator(
                    name=self.name + " comparator", suppliers=suppliers
                )
                self.supplier.is_ghost_source = True
            else:
                self.supplier = suppliers[0]
        else:
            self.supplier = suppliers
            suppliers = [suppliers]
        self.suppliers = suppliers

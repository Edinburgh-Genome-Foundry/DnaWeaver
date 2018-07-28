from ..DnaQuote import DnaQuote
from ..SequenceDecomposer import SequenceDecomposer, NoSolutionFoundError
# from .shortest_path_algorithms import NoSolutionFoundError
from .DnaSource import DnaSource
from .DnaSourcesComparator import DnaSourcesComparator
from ..DnaAssemblyMethod import (BuildAGenomeAssemblyMethod,
                                 GibsonAssemblyMethod,
                                 GoldenGateAssemblyMethod)
from ..OverhangSelector import TmOverhangSelector, FixedSizeOverhangSelector

class DnaAssemblyStation(DnaSource):
    """DNA Assembly stations assemble together DNA fragments using a specific
    assembly method.

    Parameters
    ----------

    name
      Name of the station (appears on reports)

    assembly_method
      AnDnaAssemblyMethod object specifying how the fragments are assembled,
      what sequences can be assembled, what fragments can be used, etc.

    dna_source

    """
    class_description = "DNA assembly station"
    operation_type = "assembly"
    report_fa_symbol = u""
    report_fa_symbol_plain = "flask"
    report_color = "#eeeeff"

    def __init__(self, name, assembly_method, dna_source, memoize=False,
                 decomposer_class=None, a_star_auto_multiplier=2,
                 **solve_kwargs):
        self.name = name
        self.assembly_method = assembly_method
        self.set_suppliers(dna_source)
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
        self.min_basepair_price = self.dna_source.min_basepair_price
        if solve_kwargs.get('a_star_factor', None) == 'auto':
            solve_kwargs['a_star_factor'] = 2 * self.min_basepair_price
        self.solve_kwargs = solve_kwargs

    def get_quote_for_sequence_segment(self, sequence, segment,
                                       max_lead_time=None, **kwargs):
        """Return the cost of the segment

        Is used as the "cost" function for a segment during decomposition
        optimization

        """
        fragment_to_order = self.assembly_method.compute_sequence_fragment(
            sequence, segment, **kwargs
        )
        return self.dna_source.get_quote(fragment_to_order,
                                         max_lead_time=max_lead_time)

    def get_assembly_plan_from_cuts(self, sequence, cuts, max_lead_time=None):
        """Return a plan {(seg,ment): quote, ...} based on the cut positions.
        """
        cuts = sorted(cuts)
        return {
            segment: self.get_quote_for_sequence_segment(
                sequence, segment,
                max_lead_time=max_lead_time
            )
            for segment in zip(cuts, cuts[1:])
        }

    def get_assembly_plan_for_sequence(self, sequence, max_lead_time=None,
                                       coarse_grain=None, fine_grain=None,
                                       a_star_factor=0, logger=None):
        """Return the plan {(seg,ment): quote, ...} of the optimal strategy
        for the sequence's decomposition."""
        def segment_score(segment):
            quote = self.get_quote_for_sequence_segment(
                sequence, segment, max_lead_time=max_lead_time
            )
            if not quote.accepted:
                return -1
            else:
                return quote.price
        assembly = self.assembly_method
        if coarse_grain is None:
            coarse_grain = self.solve_kwargs.get("coarse_grain", 1)
        if fine_grain is None:
            fine_grain = self.solve_kwargs.get("fine_grain", 1)

        self.decomposer = self.decomposer_class(
            sequence_length=len(sequence),
            segment_score_function=segment_score,
            cut_location_constraints=[
                cs(sequence)
                for cs in assembly.cut_location_constraints
            ],
            segment_constraints=[
                cs(sequence)
                for cs in assembly.segment_constraints
            ],
            min_segment_length=assembly.min_segment_length,
            max_segment_length=assembly.max_segment_length,
            forced_cuts=assembly.force_cuts(sequence),
            suggested_cuts=assembly.suggest_cuts(sequence),
            coarse_grain=coarse_grain,
            fine_grain=fine_grain,
            logger=logger,
            bar_prefix="[%s] " % self.name,
            a_star_factor=a_star_factor,
            path_size_limit=assembly.max_fragments,
            cuts_set_constraints=[
                cs(sequence)
                for cs in assembly.cuts_set_constraints
            ]
        )
        best_cuts = self.decomposer.compute_optimal_cuts()
        return self.get_assembly_plan_from_cuts(sequence, best_cuts,
                                                max_lead_time=max_lead_time)

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

        if (max_lead_time is not None):
            max_lead_time = max_lead_time - self.extra_time
            if max_lead_time < 0:
                return DnaQuote(self, sequence, accepted=False,
                                message="Lead time limit too short")
        try:
            assembly_plan = self.get_assembly_plan_for_sequence(
                sequence, max_lead_time=max_lead_time, **self.solve_kwargs
            )
        except NoSolutionFoundError:
            return DnaQuote(self, sequence, accepted=False,
                            message="No solution found !")

        # A solution has been found ! Now compute overall time and lead time.
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
            "dna source": self.dna_source.name,
            "solver parameters": self.solve_kwargs,
        }
        result.update({
            ("assembly method %s" % k): v
            for k, v in self.assembly_method.dict_description().items()
        })
        return result

    @staticmethod
    def from_dict(data):
        if data['method'] == 'golden_gate_assembly':
            gc_range = data.get('overhang_gc_range', [0, 1])
            method = GoldenGateAssemblyMethod(
                min_overhangs_gc=gc_range[0],
                max_overhangs_gc=gc_range[1],
                enzyme=data['enzyme']
            )
        elif data['method'] in ['gibson_assembly', 'yeast_recombination']:
            if data['overhang_type'] == 'tm':
                min_oh_size, max_oh_size = data['overhang_size_range']
                min_tm, max_tm = data['tm_range']
                overhang_selector = TmOverhangSelector(
                    min_size=min_oh_size, max_size=max_oh_size,
                    min_tm=min_tm, max_tm=max_tm
                )
            else:
                overhang_selector = FixedSizeOverhangSelector(
                    overhang_size=data['overhang_size']
                )

            method = GibsonAssemblyMethod(
                overhang_selector=overhang_selector
            )

        return DnaAssemblyStation(
            name=data['name'],
            dna_source=data['suppliers'],
            assembly_method=method,
            coarse_grain=10,
            fine_grain=2,
            logger='bars'
        )

    def set_suppliers(self, suppliers):
        if hasattr(suppliers, '__iter__'):
            if len(suppliers) > 1:
                self.dna_source = DnaSourcesComparator(
                    name=self.name + ' comparator', suppliers=suppliers)
            else:
                self.dna_source = suppliers[0]
        else:
            self.dna_source = suppliers

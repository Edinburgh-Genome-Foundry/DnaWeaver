from .biotools import (reverse_complement, find_enzyme_sites,
                       gc_content_to_tm)
from .OverhangSelector import TmOverhangSelector
from .tools import memoize
import itertools


class AssemblyMethod:
    """General class for assembly methods.

    All assembly methods can store the following attributes:

    Parameters
    ----------

    duration
      Duration required for the assembly (e.g. an estimated upper bound). The
      time unit is left at the choice of the user.

    cut_location_constraints
      List or tuple of functions `(sequence, int) -> bool` which return for
      a sequence and an index (cut location) whether the location should be
      considered as a cutting site compatible with this assembly method (True)
      or not.
      The locations considered in fine are the locations which pass every
      constraint in the `cut_location_contraints` list.

    segment_constraints
      List or tuple of functions `(sequence, start, end) -> bool` which returns
      for a subsegment `(start, end)` whether the segment is a valid segment
      for the assembly.
      The segments considered in fine are the segments which pass every
      constraint in the `segments_constraints` list.

    min_segment_length
      Minimal length of the fragments that this assembly method allows

    max_segment_length
      Maximal length of the fragments that this assembly method allows


    force_cuts
      A function `sequence->[i1, i2, i3...]` which for a given sequence returns
      forced cuts locations.

    max_fragments
       Maximal number of fragments which can be assembled with this method.

    sequence_constraints
      General constraints on the sequence so that it can be built with this
      assembly method. It is a list of functions `seq->bool`. If one of
      these constraints returns False the sequence is refused.

    reference
      Reference to e.g. a paper or a protocol describing the method.
    """
    name = "None"

    def __init__(self, duration=0, cost=0, reference=None,
                 cut_location_constraints=(),
                 segment_constraints=(), min_segment_length=None,
                 max_segment_length=None, force_cuts=(), suggest_cuts=(),
                 max_fragments=None, sequence_constraints=(),
                 cuts_set_constraints=()):
        self.duration = duration
        self.cost = cost
        self.cut_location_constraints = list(cut_location_constraints)
        self.segment_constraints = list(segment_constraints)
        self.min_segment_length = min_segment_length
        self.max_segment_length = max_segment_length
        self.sequence_constraints = list(sequence_constraints)
        self.max_fragments = max_fragments
        self.cuts_set_constraints = list(cuts_set_constraints)
        self.reference = reference

        if callable(suggest_cuts):
            self.suggest_cuts = suggest_cuts
        else:
            # means the cuts are a constant
            self.suggest_cuts = lambda *a: suggest_cuts

        if callable(force_cuts):
            self.force_cuts = force_cuts
        else:
            # means the forced cuts are a constant
            self.force_cuts = lambda *a: force_cuts

    def dict_description(self):
        result = {
            "name": self.name,
            "minimum segment length": self.min_segment_length,
            "maximum segment length": self.max_segment_length,
            "maximal number of fragments": self.max_fragments,
            "reference": self.reference,
            "duration": self.duration,
            "cost": self.cost
        }
        result.update(self.additional_dict_description())
        return result

    def additional_dict_description(self):
        return {}


class OverlapingAssemblyMethod(AssemblyMethod):
    """General class for all overlaping assembly methods.

    Parameters
    ----------

    homology_arm_length
      Length of the homology arm, or "overhang". A length of L means that
      consecutive segments will overlap by 2*L

    """
    name = "Overlaping Assembly"

    def __init__(self, overhang_selector, **properties):
        AssemblyMethod.__init__(self, **properties)
        self.overhang_selector = overhang_selector
        self.cut_location_constraints.append(
            overhang_selector.location_filter_method)
        self.compute_sequence_fragment = overhang_selector.compute_sequence_fragment


    # def compute_fragment_sequence(self, sequence):
    #     """Return the segment's sequence with flanking sequences.
    #
    #     Parameters
    #     ----------
    #
    #     segment
    #       A pair of integers (start, end) delimiting the subfragment
    #       sequence[start:stop]
    #
    #     sequence
    #       An "ATGC" DNA sequence string
    #
    #     """
    #
    #
    #     def f(segment):
    #
    #     return self.overhang_selector.compute_fragment_sequence(sequence,
    #                                                             segment)
    #     L = len(sequence)
    #     start, end = segment
    #     return sequence[max(0, start - self.homology_arm_length):
    #                     min(L, end + self.homology_arm_length)]


class GibsonAssemblyMethod(OverlapingAssemblyMethod):
    """Gibson Assembly Method. Just another overlap-method"""
    name = "Gibson Assembly"




class BuildAGenomeAssemblyMethod(OverlapingAssemblyMethod):
    """The Build-a-Genome Assembly Method. Just another overlap-method"""
    name = "Build-a-Genome"


class GoldenGateAssemblyMethod(OverlapingAssemblyMethod):
    """The Golden Gate Assembly Method.

    This method adds overhangs with a Type IIS REase site to the segments.

    Parameters
    ----------

    left_overhang
      An ATGC DNA sequence representing the left overhang, in 5' -3'.
      For practicality, this can contain "[BsaI]", "[BsmBI]", "[BbsI]", which
      will be replaced by their restriction sites


    right_overhang
      An ATGC DNA sequence representing the right overhang, in 5'-3'.
      Be careful, I said 5'-3' !!!
      If left to None, will be equal to left_overhang.

    """
    name = "Golden Gate Assembly"

    enzymes_dict = {
        "BsaI": "GGTCTC",
        "BsmBI": "CGTCTC",
        "BbsI": "GAAGAC",
    }

    def __init__(self, enzyme="BsaI", wildcard_basepair="A",  left_overhang="",
                 right_overhang="", refuse_sequences_with_enzyme_site=True,
                 min_overhangs_gc=0, max_overhangs_gc=1,
                 min_overhangs_differences=1,  **props):
        if enzyme not in self.enzymes_dict:
            return ValueError("Enzyme should be one of %s" %
                              self.enzymes_dict.keys())

        # def location_has_valid_gc_content(sequence):
        #     def f(cut_location):
        #         gc = gc_content(get_overhang(sequence, cut_location))
        #         return (self.min_overhangs_gc < gc < self.max_overhangs_gc)
        #     return f
        # self.cut_location_constraints.append(location_has_valid_gc_content)
        self.min_overhangs_gc = min_gc = min_overhangs_gc
        self.max_overhangs_gc = max_gc = max_overhangs_gc
        overhang_selector = TmOverhangSelector(
            min_size=4, max_size=4,
            min_tm=gc_content_to_tm(4, min_gc),
            max_tm=gc_content_to_tm(4, max_gc)
        )
        OverlapingAssemblyMethod.__init__(self, overhang_selector, **props)
        self.enzyme = enzyme
        self.enzyme_site = self.enzymes_dict[enzyme]
        enzyme_site_plus_basepair = self.enzyme_site + wildcard_basepair
        self.left_overhang = left_overhang + enzyme_site_plus_basepair
        self.right_overhang = right_overhang + enzyme_site_plus_basepair
        self.right_overhang_rev = reverse_complement(self.right_overhang)
        self.min_overhangs_differences = min_overhangs_differences
        self.refuse_sequences_with_enzyme_site = \
            refuse_sequences_with_enzyme_site

        # CUTS LOCATION CONSTRAINT BASED ON GC CONTENT

        if refuse_sequences_with_enzyme_site:
            def no_site_in_sequence(sequence):
                sites = find_enzyme_sites(sequence, enzyme_name=self.enzyme)
                return sites == []
            self.sequence_constraints.append(no_site_in_sequence)

        # CUTS LOCATION CONSTRAINT BASED ON GC CONTENT



        # CUTS SET CONSTRAINT: ALL OVERHANGS MUST BE COMPATIBLE

        def overhangs_are_compatible(o1, o2):
            diffs = len([a for a, b in zip(o1, o2) if a != b])
            if diffs >= self.min_overhangs_differences:
                rev_o2 = reverse_complement(o2)
                rev_diffs = len([a for a, b in zip(o1, rev_o2) if a != b])
                return rev_diffs >= self.min_overhangs_differences
            return False

        overhangs_are_compatible = memoize(overhangs_are_compatible)

        def all_overhangs_are_compatible(sequence):
            def f(cut_locations):
                overhangs = sorted([
                    overhang_selector.compute_overhang_sequence(sequence,
                                                                cut_location)
                    for cut_location in cut_locations
                ])
                return all([
                    overhangs_are_compatible(o1, o2)
                    for o1, o2 in itertools.combinations(overhangs, 2)
                ])
            return f

        self.cuts_set_constraints.append(all_overhangs_are_compatible)



    # def compute_fragment_sequence(self, sequence, segment):
    #     """Return the segment's sequence with flanking sequences for
    #
    #     Parameters
    #     ----------
    #
    #     segment
    #       A pair of integers (start, end) delimiting the subfragment
    #       sequence[start:stop]
    #
    #     sequence
    #       An "ATGC" DNA sequence string
    #
    #     """
    #     L = len(sequence)
    #     start, end = segment
    #     segment = sequence[max(0, start - 2): min(L, end + 2)]
    #     return self.left_overhang + segment + self.right_overhang_rev

    def additional_dict_description(self):
        return {
            "enzyme": self.enzyme,
            "left overhang": self.left_overhang,
            "right overhang": self.right_overhang,
            "refuse sequences with enzyme site":
                str(self.refuse_sequences_with_enzyme_site),
            "overhangs gc content": "%d-%d%%" % (100*self.min_overhangs_gc,
                                                 100*self.max_overhangs_gc),
            "overhangs differences": self.min_overhangs_differences
        }

import itertools
from ..biotools import (
    reverse_complement,
    gc_content_to_tm,
    find_enzyme_sites,
    get_sequence_topology,
)

from ..tools import memoize
from ..SegmentSelector import TmSegmentSelector
from .OverlapingAssemblyMethod import OverlapingAssemblyMethod


class GoldenGateAssemblyMethod(OverlapingAssemblyMethod):
    """The Golden Gate Assembly Method.

    This method adds overhangs with a Type IIS REase site to the segments.

    Parameters
    ----------

    left_overhang
      An ATGC DNA sequence representing the left overhang, in 5' -3'.
      For practicality, this can contain "[BsaI]", "[BsmBI]", "[BbsI]", which
      will be replaced by their restriction sites.

    right_overhang
      An ATGC DNA sequence representing the right overhang, in 5'-3'.
      Be careful, it has to be 5'-3' !!!
      If left to None, will be equal to left_overhang.
    """

    name = "Golden Gate Assembly"

    enzymes_dict = {
        "BsaI": "GGTCTC",
        "BsmBI": "CGTCTC",
        "BbsI": "GAAGAC",
        "SapI": "GCTCTTC",
    }

    def __init__(
        self,
        enzyme="BsaI",
        wildcard_basepair="A",
        left_addition="",
        right_addition="",
        refuse_sequences_with_enzyme_site=True,
        min_overhangs_gc=0,
        max_overhangs_gc=1,
        min_overhangs_differences=1,
        **props
    ):
        if enzyme not in self.enzymes_dict:
            raise ValueError("Enzyme should be one of %s" % self.enzymes_dict.keys())

        self.min_overhangs_gc = min_gc = min_overhangs_gc
        self.max_overhangs_gc = max_gc = max_overhangs_gc
        self.min_overhangs_differences = min_overhangs_differences

        self.enzyme = enzyme
        self.enzyme_site = self.enzymes_dict[enzyme]
        enzyme_site_plus_basepair = self.enzyme_site + wildcard_basepair
        self.left_addition = left_addition + enzyme_site_plus_basepair
        self.right_addition = (
            reverse_complement(enzyme_site_plus_basepair) + right_addition
        )
        self.refuse_sequences_with_enzyme_site = refuse_sequences_with_enzyme_site
        self.overhang_size = 3 if self.enzyme == "SapI" else 4

        overhang_selector = TmSegmentSelector(
            min_size=self.overhang_size,
            max_size=self.overhang_size,
            min_tm=gc_content_to_tm(self.overhang_size, min_gc),
            max_tm=gc_content_to_tm(self.overhang_size, max_gc),
            left_addition=self.left_addition,
            right_addition=self.right_addition,
        )
        OverlapingAssemblyMethod.__init__(self, overhang_selector, **props)

        # CUTS LOCATION CONSTRAINT BASED ON GC CONTENT

        if refuse_sequences_with_enzyme_site:

            def no_site_in_sequence(sequence):
                sites = find_enzyme_sites(sequence, enzyme_name=self.enzyme)
                return sites == []

            self.sequence_constraints.append(no_site_in_sequence)

        # DO NOT CUT AT PALINDROMIC REGIONS

        def no_cut_at_palyndromic_locations(sequence):
            def no_palyndrom_filter(i):
                s = overhang_selector.compute_segment_around_index(sequence, i)
                rev_s = reverse_complement(s)
                rev_diffs = len([a for a, b in zip(s, rev_s) if a != b])
                assert len(s) == self.overhang_size
                return rev_diffs >= self.min_overhangs_differences

            return no_palyndrom_filter

        self.cut_location_constraints.append(no_cut_at_palyndromic_locations)

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
            topology = get_sequence_topology(sequence, "linear")

            def constraint(cut_locations):
                cut_overhangs = {
                    cut_location: overhang_selector.compute_segment_around_index(
                        sequence, cut_location
                    )
                    for cut_location in cut_locations
                }
                cut_pairs = list(itertools.combinations(cut_locations, 2))
                if topology == "circular":
                    cut_pairs.remove((0, len(sequence)))

                return all(
                    [
                        overhangs_are_compatible(cut_overhangs[c1], cut_overhangs[c2])
                        for c1, c2 in cut_pairs
                    ]
                )

                # overhangs = sorted(
                #     [
                #         overhang_selector.compute_segment_around_index(
                #             sequence, cut_location
                #         )
                #         for cut_location in cut_locations
                #     ]
                # )
                # return all(
                #     [
                #         overhangs_are_compatible(o1, o2)
                #         for o1, o2 in itertools.combinations(overhangs, 2)
                #     ]
                # )

            return constraint

        self.cuts_set_constraints.append(all_overhangs_are_compatible)

    def additional_dict_description(self):
        return {
            "enzyme": self.enzyme,
            "left addition": self.left_addition,
            "right addition": self.right_addition,
            "refuse sequences with enzyme site": str(
                self.refuse_sequences_with_enzyme_site
            ),
            "overhangs gc content": "%d-%d%%"
            % (100 * self.min_overhangs_gc, 100 * self.max_overhangs_gc),
            "overhangs differences": self.min_overhangs_differences,
        }

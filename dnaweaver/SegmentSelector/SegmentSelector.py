""" This module contains overhangs selectors for methods such as
Gibson Assembly, Golden Gate assembly, recombination in yeast."""

from ..biotools import get_sequence_topology
import numpy as np


class SegmentSelector:
    """Base class for segment selectors such as TmSegmentSelector.

    These selectors return an ideal segment subsegment around a given sequence
    location. They can also be used to filter out locations in the sequence as
    non-potential cutting sites for sequence decomposition.
    """

    def location_filter_method(self, sequence):
        """Return a filter function f(location) => True/False.

        The result is True if the location is a valid cutting site."""
        if self.has_location_filter:

            def f(index):
                return self.filter_location(sequence, index)

            return f
        else:
            return None

    @property
    def has_location_filter(self):
        return hasattr(self, "filter_location")

    def compute_fragment_for_sequence_segment(self, sequence, segment):
        """Compute the fragment as (sequence region + flank segments)."""
        start, end = segment
        topology = get_sequence_topology(sequence)
        if topology == "circular":
            half_homology = int(np.ceil(self.max_homology_size / 2))
            if start <= half_homology:
                return self.compute_fragment_for_sequence_segment(
                    sequence=sequence[-half_homology:] + sequence,
                    segment=(start + half_homology, end + half_homology),
                )
            if end >= len(sequence) - half_homology:
                return self.compute_fragment_for_sequence_segment(
                    sequence=sequence + sequence[:half_homology], segment=segment,
                )
        if start == 0:
            fragment_start = 0
        else:
            loc = self.compute_segment_location(sequence, start)
            if loc is None:
                subseq = sequence[:20] + "..."
                raise ValueError(
                    "loc was None around start %s for sequence %s (%dbp)"
                    % (end, subseq, len(sequence))
                )
            fragment_start = loc[0]
        if end == len(sequence):
            fragment_end = len(sequence)
        else:
            loc = self.compute_segment_location(sequence, end)
            if loc is None:
                subseq = sequence[:20] + "..."
                raise ValueError(
                    "loc was None around end %s for sequence %s (%dbp)"
                    % (end, subseq, len(sequence))
                )
            fragment_end = loc[1]
        fragment = sequence[fragment_start:fragment_end]
        return self.left_addition + fragment + self.right_addition

    def compute_segment_around_index(self, sequence, index, topology="linear"):
        """Return the sequence of the selected segment at the given index."""

        if get_sequence_topology(sequence) == "circular":
            half_homology = int(np.ceil(self.max_homology_size / 2))
            if index <= half_homology:
                return self.compute_segment_around_index(
                    sequence=sequence[-half_homology:] + sequence,
                    index=index + half_homology,
                )
            if index >= len(sequence) - half_homology:
                return self.compute_segment_around_index(
                    sequence=sequence + sequence[:half_homology], index=index,
                )
        start, end = self.compute_segment_location(sequence, index)
        return sequence[start:end]

    @staticmethod
    def get_segment_coordinates(center, segment_length, sequence_length):
        """Return max(0, c - s/2) - min(L, c + L/2).

        Where c=center, s=segment_length, L=sequence_length.
        """
        half = int(segment_length / 2)
        start = max(0, min(center - half, sequence_length - segment_length))
        end = start + segment_length
        return start, end

    def __repr__(self):
        return str(self)

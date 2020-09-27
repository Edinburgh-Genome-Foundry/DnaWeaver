import re
from io import StringIO
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO
from .SequenceString import SequenceString
import numpy as np


def get_sequence_topology(sequence, default="linear"):
    if isinstance(sequence, SequenceString):
        return sequence.metadata.get("topology", default)
    else:
        return default


def gc_content(sequence, window_size=None):
    """Compute global or local GC content.

    Parameters
    ----------

    sequence
      An ATGC DNA sequence (upper case!).

    window_size
      If provided, the local GC content for the different sliding windows of
      this size is returned, else the global GC content is returned.

    Returns
    --------

      A number between 0 and 1, indicating the proportion
      of GC content. If window_size is provided, returns
      a list of len(sequence)-window_size values indicating
      the local GC contents (sliding-window method). The i-th value
      indicates the GC content in the window [i, i+window_size].
    """
    # The code is a little cryptic but speed gain is 300x
    # compared with pure-python string operations
    arr = np.frombuffer((sequence + "").encode(), dtype="uint8")
    arr_GCs = (arr == 71) | (arr == 67)  # 67=C, 71=G

    if window_size is None:
        return 1.0 * arr_GCs.sum() / len(sequence)
    else:
        cs = np.cumsum(arr_GCs)
        a = cs[window_size - 1 :]
        b = np.hstack([[0], cs[:-window_size]])
        return 1.0 * (a - b) / window_size


def find_enzyme_sites(sequence, enzyme_name, padding=0, padding_nuc="A"):
    padding = padding * padding_nuc
    if hasattr(sequence, "seq"):
        sequence = str(sequence.seq)
    sequence = str(sequence)
    sequence = Seq(padding + sequence + padding)
    topology = get_sequence_topology(sequence, default="linear")
    linear = topology == "linear"
    return Restriction.__dict__[enzyme_name].search(sequence, linear=linear)


def gc_content_to_tm(seq_length, gc_content, gc_tm=4, nongc_tm=2):
    return seq_length * (gc_content * gc_tm + (1 - gc_content) * nongc_tm)

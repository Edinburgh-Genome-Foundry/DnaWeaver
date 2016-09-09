from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import numpy as np
import subprocess
import tempfile
import time
import os
import re


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.
    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.
    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def random_dna_sequence(length, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def blast_sequence(sequence, blast_db=None, subject=None, word_size=4,
                   perc_identity=80, num_alignments=1000, num_threads=3,
                   use_megablast=True, ungapped=True):
    """Return a Biopython BLAST record of the given sequence BLASTed
    against the provided database.

    Parameters
    ----------

    sequence
      An ATGC sequence

    Examples
    --------

    >>> blast_record = blast_sequence("ATTGTGCGTGTGTGCGT", "blastdb/ecoli")
    >>> for alignment in blast_record.alignments:
    >>>     for hit in alignment.hsps:
    >>>         print (hit.identities)
    """

    xml_file, xml_name = tempfile.mkstemp(".xml")
    fasta_file, fasta_name = tempfile.mkstemp(".fa")
    with open(fasta_name, "w+") as f:
        f.write(">seq\n" + sequence)


    if subject is not None:
        close_subject = True
        if not subject.endswith(".fa"):
            remove_subject = True
            subject_file, fasta_subject_name = tempfile.mkstemp(".fa")
            with open(fasta_subject_name, "w+") as f:
                f.write(">subject\n" + subject)
            subject = fasta_subject_name
        else:
            remove_subject = False
    else:
        close_subject = False
    print (fasta_name, subject)

    p = subprocess.Popen([
        "blastn", "-out", xml_name,
        "-outfmt", "5",
        "-num_alignments", str(num_alignments),
        "-query", fasta_name] +
        (["-db", blast_db] if blast_db is not None
         else ['-subject', subject]) +
        (["-ungapped"] if ungapped else []) +
        (["-task", "megablast"] if use_megablast else []) + [
        "-word_size", str(word_size),
        "-num_threads", str(num_threads),
        "-perc_identity", str(perc_identity)
    ], close_fds=True, stderr=subprocess.PIPE)
    res, blast_err = p.communicate()
    p.wait()
    print(("wwwww", xml_name))
    error = None
    for i in range(3):
        try:
            with open(xml_name, "r") as f:
                res = list(NCBIXML.parse(f))
                os.fdopen(xml_file, 'w').close()
                os.fdopen(fasta_file, 'w').close()
                os.remove(xml_name)
                os.remove(fasta_name)

                if close_subject:
                    open(subject, 'w').close()
                    if remove_subject:
                        os.remove(subject)
                if len(res) == 1:
                    return res[0]
                else:
                    return res
            break
        except ValueError as err:
            error = err
            time.sleep(0.1)
    else:
        raise ValueError("Problem reading the blast record: " + str(error))




def largest_common_substring(query, target, max_overhang):
    """Return the largest common substring between `query` and `target`.

    If the common substring is too much smaller than `query` False is returned,
    else the location `(start, end)` of the substring in `target` is returned.

    Parameters:
    -----------

    query (str)
      The sequence to be found in target (minus some overhangs possibly)

    target (str)
      The sequence in which to find `query`

    max_overhang
      Maximal size allowed for the flanking regions of `query` that would
      not be contained in `target`.

    Notes:
    ------

    This is intended for finding whether `query` can be extracted from `target`
    using PCR. See the PcrOutStation implementation in DNASource.py.

    """
    start, end = max_overhang, len(query) - max_overhang
    if query[start:end] not in target:
        return False
    while (start >= 0) and (query[start:end] in target):
        start -= 1
    start += 1
    while (end < len(query)) and (query[start:end] in target):
        end += 1
    end -= 1
    return start, end


def gc_content(sequence, window_size=None):
    """Compute global or local GC content.

    Parameters
    ----------

    sequence
      An ATGC DNA sequence (upper case!)

    window_size
      If provided, the local GC content for the different sliding windows of
      this size is returned, else the global GC content is returned.

    Returns
    --------

      A number between 0 and 1 indication the proportion
      of GC content. If window_size is provided, returns
      a list of len(sequence)-window_size values indicating
      the local GC contents (sliding-window method). The i-th value
      indicates the GC content in the window [i, i+window_size]
    """
    # The code is a little cryptic but speed gain is 300x
    # compared with pure-python string operations

    arr = np.fromstring(sequence + "", dtype="uint8")
    arr_GCs = (arr == 71) | (arr == 67)  # 67=C, 71=G

    if window_size is None:
        return 1.0 * arr_GCs.sum() / len(sequence)
    else:
        cs = np.cumsum(arr_GCs)
        a = cs[window_size - 1:]
        b = np.hstack([[0], cs[:-window_size]])
        return 1.0 * (a - b) / window_size


def no_pattern_constraint(pattern, is_regex=False, with_revcomp=True):
    if is_regex:
        cm_pattern = re.compile(pattern)

        def not_in_strand(sequence):
            return (cm_pattern.search(sequence) is not None)
    else:
        def not_in_strand(sequence):
            return pattern not in sequence
    if with_revcomp:
        def not_in_both_strands(sequence):
            if not_in_strand(sequence):
                if not_in_strand(reverse_complement(sequence)):
                    return True
            return False
        return not_in_both_strands
    else:
        return not_in_strand

import tempfile
import time
import os
import subprocess
from Bio.Blast import NCBIXML
import numpy as np
from .sequence_operations import sequence_to_atgc


def blast_sequence(
    sequence,
    blast_db=None,
    subject=None,
    word_size=4,
    perc_identity=80,
    num_alignments=1000,
    num_threads=3,
    use_megablast=True,
    ungapped=True,
):
    """Return a Biopython BLAST record of the given sequence BLASTed
    against the provided database.

    Parameters
    ----------

    sequence
      An ATGC sequence.

    blast_db
      Path to a BLAST database.

    subject
      Either a path to a fasta (.fa) file or an ATGC string. Subject to blast
      against.

    word_size
      Word size to use in the blast.

    perc_identity
      Minimal percentage of identical nucleotides in a match for it to be kept.

    num_alignments
      Number of alignments to keep.

    num_threads
      Number of threads for the BLAST.

    use_megablast
      Whether to use Megablast.

    ungapped
      No-gaps matches only ?

    Examples
    --------

    >>> blast_record = blast_sequence("ATTGTGCGTGTGTGCGT", "blastdb/ecoli")
    >>> for alignment in blast_record.alignments:
    >>>     for hit in alignment.hsps:
    >>>         print (hit.identities)
    """

    xml_file, xml_name = tempfile.mkstemp(".xml")
    fasta_file, fasta_name = tempfile.mkstemp(".fa")
    sequence = sequence_to_atgc(sequence)
    with open(fasta_name, "w+") as f:
        f.write(">seq\n" + sequence)

    if subject is not None:
        close_subject = True
        if not subject.endswith(".fa"):
            remove_subject = True
            _subject_file, fasta_subject_name = tempfile.mkstemp(".fa")
            with open(fasta_subject_name, "w+") as f:
                f.write(">subject\n" + subject)
            subject = fasta_subject_name
        else:
            remove_subject = False
    else:
        close_subject = False
    p = subprocess.Popen(
        [
            "blastn",
            "-out",
            xml_name,
            "-outfmt",
            "5",
            "-num_alignments",
            str(num_alignments),
            "-query",
            fasta_name,
        ]
        + (["-db", blast_db] if blast_db is not None else ["-subject", subject])
        + (["-ungapped"] if ungapped else [])
        + (["-task", "megablast"] if use_megablast else [])
        + [
            "-word_size",
            str(word_size),
            "-num_threads",
            str(num_threads),
            "-dust",
            "no",
            "-evalue",
            "0.01",
            "-perc_identity",
            str(perc_identity),
        ],
        close_fds=True,
        stderr=subprocess.PIPE,
    )
    res, _blast_err = p.communicate()
    p.wait()
    error = None
    for i in range(3):
        try:
            with open(xml_name, "r") as f:
                res = list(NCBIXML.parse(f))
                os.fdopen(xml_file, "w").close()
                os.fdopen(fasta_file, "w").close()
                os.remove(xml_name)
                os.remove(fasta_name)

                if close_subject:
                    open(subject, "w").close()
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


def make_blast_db(fasta_input, target):
    proc = subprocess.Popen(
        ["makeblastdb", "-in", fasta_input, "-dbtype", "nucl", "-out", target]
    )
    proc.wait()


def perfect_match_locations_in_hsp(hsp, span_cutoff=10):
    """Return the locations of perfect matches in a BLAST HSP.

    Only locations with a span above span_cutoff are kept.
    """
    if hsp.align_length < span_cutoff:
        return []
    arr = np.frombuffer(hsp.match.encode(), dtype="uint8")
    indices = [0] + list((arr != 124).nonzero()[0]) + [len(arr)]
    return [
        (start + hsp.query_start, end + hsp.query_start)
        for start, end in zip(indices, indices[1:])
        if end - start >= span_cutoff
    ]


def largest_common_substring(query, target, max_overhang):
    """Return the largest common substring between `query` and `target`.

    Find the longest substring of query that is contained in target.

    If the common substring is too much smaller than `query` False is returned,
    else the location `(start, end)` of the substring in `target` is returned.

    Parameters:
    -----------

    query (str)
      The sequence to be found in target (minus some overhangs possibly).

    target (str)
      The sequence in which to find `query`.

    max_overhang
      Maximal size allowed for the flanking regions of `query` that would
      not be contained in `target`.

    Examples
    --------

    >>> seqA = '-----oooooooo'
    >>> seqB = 'oooooo-----tttt'
    >>> largest_common_substring(seqA, seqA, 80) # == (0, 12)
    >>> largest_common_substring(seqA, seqB, 80) # == (5, 11)

    Notes:
    ------

    This is intended for finding whether `query` can be extracted from `target`
    using PCR. See the PcrExtractionStation implementation in DnaSupplier.py.
    """
    # The trick here is to start with the central region of "query".
    # This region is initially as small as max_overhang allows, and it is
    # progressively expanded on the sides
    max_overhang = min(max_overhang, int(len(query) / 2))
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

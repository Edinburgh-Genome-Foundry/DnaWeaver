from Bio import NCBIXML
from Bio.Seq import Seq
import subprocess


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


def blast_sequence(sequence, blast_db, word_size=4, perc_identity=80,
                   num_alignments=1000, num_threads=3):
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

    p = subprocess.Popen([
        "blastn", "-out", xml_name,
        "-outfmt", "5",
        "-num_alignments", str(num_alignments),
        "-query", fasta_name,
        "-db", blast_db,
        "-word_size", str(word_size),
        "-num_threads", str(num_threads),
        "-perc_identity", str(perc_identity)
    ], close_fds=True)
    p.communicate()
    p.wait()
    for i in range(3):
        try:
            with open(xml_name, "r") as f:
                blast_record = NCBIXML.read(f)
            break
        except ValueError:
            time.sleep(0.1)
    else:
        raise ValueError("Problem reading the blast record.")

    os.fdopen(xml_file, 'w').close()
    os.fdopen(fasta_file, 'w').close()

    return blast_record

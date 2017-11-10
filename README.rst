DnaWeaver
==========

(documentation in progress, come back later !)

DnaWeaver is a Python library to find optimal strategies to assemble large
synthetic DNA fragments.

A DNA assembly problem is defined by the sequence to be assembled, and a supply
network of DNA sources (such as commercial offers, parts libraries, or assembly
stations) which produce or assemble sub-fragments of DNA.

Given such a problem, DnaWeaver produces [a report](example) describing an
optimized assembly plan: what sub-fragments should be ordered, what sub-fragments
should be obtained from an existing construct/genome, what cloning methods
should be used for each assembly step, etc.

DnaWeaver was written with versatility and extensibility in mind:
each DNA source and assembly method can be customized, and assembly plans can
be optimized with respect to total price, overall duration of the assembly,
or assembly success probabilities.

DnaWeaver can also export the result as interactive widgets for web applications, and
as JSON for automated assembly platforms.

Example of use
---------------

In the following example we compute an assembly plan for a 10000bp sequence,
where the sub-fragments can be ordered from two companies:

- **CheapDNA** produces fragments for 10c/bp, at the condition that they do not
  contain any BsaI site (GGTCTC) and that they are under 4000bp in size.
- **DeluxeDNA** produces any fragment under 3000bp for 20c/bp.

Fragments are assembled using Gibson Assembly with 40bp overlap between segments.

Here is the Python code to solve the problem with DnaWeaver:

.. code:: python

    from dnaweaver import *

    cheap_dna_offer = CommercialDnaOffer(
        name="CheapDNA.com",
        sequence_constraints=[
            no_pattern_constraint("GGTCTC"),
            lambda seq: len(seq) < 4000
        ],
        pricing=lambda seq: 0.10 * len(seq),
    )

    deluxe_dna_offer = CommercialDnaOffer(
        name="DeluxeDNA.com",
        sequence_constraints=[lambda seq: len(seq) < 3000],
        pricing=(lambda seq: 0.20 * len(seq)),
    )

    assembly_station = DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=GibsonAssemblyMethod(homology_arm_length=20,
                                             min_segment_length=500,
                                             max_segment_length=4000),
        dna_source=DnaSourcesComparator([cheap_dna_offer, deluxe_dna_offer]),
        nucleotide_resolution=10,
        refine_resolution=False
    )

    sequence = random_dna_sequence(10000, seed=123)
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)

    print (quote.assembly_step_summary())

Result:

.. code:: bash

    Ordering plan:
    (0, 2005): From CheapDNA.com, price 202.50, lead_time 10.0
    (2005, 4020): From CheapDNA.com, price 205.50, lead_time 10.0
    (4020, 6025): From DeluxeDNA.com, price 409.00, lead_time 5.0
    (6025, 10000): From CheapDNA.com, price 399.50, lead_time 10.0


See the examples section for more complete examples involving different sources,
multiple assembly methods, and complex biological constraints.


Installation
-------------

You can install DnaWeaver through PIP
::
    sudo pip install dnaweaver

Alternatively, you can unzip the sources in a folder and type
::
    sudo python setup.py install


License = MIT
--------------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

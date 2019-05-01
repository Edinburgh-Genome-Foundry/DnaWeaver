.. raw:: html

    <p align="center">
    <img alt="DNA Weaver Logo" title="DNA Weaver Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaWeaver/master/docs/_static/images/title.png" width="500">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaWeaver.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaWeaver
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaWeaver/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaWeaver?branch=master


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

    import dnaweaver as dw

    cheap_dna_offer = dw.CommercialDnaOffer(
        name="CheapDNA.com",
        sequence_constraints=[
            dw.NoPatternConstraint(enzyme="BsaI"),
            dw.SequenceLengthConstraint(max_length=4000)
        ],
        pricing=dw.PerBasepairPricing(0.10),
    )

    deluxe_dna_offer = dw.CommercialDnaOffer(
        name="DeluxeDNA.com",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=3000)],
        pricing=dw.PerBasepairPricing(0.20),
    )

    assembly_station = dw.DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.TmOverhangSelector(min_tm=55, max_tm=70),
            min_segment_length=500,
            max_segment_length=4000
        ),
        dna_source=[cheap_dna_offer, deluxe_dna_offer],
        logger='bar',
        coarse_grain=20,
        fine_grain=1
    )

    sequence = dw.random_dna_sequence(10000, seed=123)
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

Also install the ncbi-blast+ package. On Ubuntu:

::
    sudo apt-get install ncbi-blast+

Reports generation needs more dependencies. Install Python dependencies with

::
    sudo pip install pandas dna_features_viewer weasyprint

Install non-python dependencies as follows on Ubuntu:
::
    sudo apt-get installbuild-essential python3-dev python3-pip \
        python3-cffi libcairo2 libpango-1.0-0 libpangocairo-1.0-0 \
        libgdk-pixbuf2.0-0 libffi-dev shared-mime-info

License = MIT
--------------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

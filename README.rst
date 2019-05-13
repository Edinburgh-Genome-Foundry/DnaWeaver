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

DNA Weaver is a Python library to find optimal strategies for assembling large
DNA constructs. Given an arbitrary sequence, DNA Weaver it will select the most
adapted commercial DNA providers, cloning methods and parts repositories
(depending on your preferences), and will design all necessary assembly
fragments and assembly steps. Try it online via the
[DNA Weaver web app](https://dnaweaver.genomefoundry.org)!

DNA Weaver was written with versatility and extensibility in mind:
each DNA source and assembly method can be customized, and assembly plans can
be optimized with respect to total price, overall duration of the assembly,
or assembly success probabilities.

How it works
------------

In DNA Weaver you first define a supply network connecting various DNA sources
(commercial providers, parts repositories, genomic DNA, and cloning stations) to
represent how DNA can be obtained in your lab or biofoundry. For instance, assume
that you routinely assemble ~10kb sequences using Gibson assembly, from fragments
obtained either commercially or from the assembly of oligonucleotides. Your
supply network then looks as follows:

.. raw:: html

    <p align="center">
    <img alt="DNA Weaver Logo" title="DNA Weaver Logo"
         src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaWeaver/master/docs/_static/images/vendor_or_oligo_assembly.png" width="250"/>
    <br /><br />
    </p>

When you submit a sequence to the main station (here, the Gibson Assembly station),
DNA Weaver will use smart sequence decomposition techniques and competitive
bidding between the different DNA sources in order to find the best possible
assembly plan.

Examples
---------

Ordering fragments from multiple vendors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following example we ask DNA Weaver for a plan to assemble a 10kb
sequence via Gibson assembly of fragments with 40bp homologies. The fragments
can be ordered from two companies: CheapDNA produces fragments for 10c/bp,
at the condition that they do not contain any BsaI site and are smaller than 3kb,
and DeluxeDNA produces produces any fragment under 4kb for 20c/bp.

We will first CheapDNA and DeluxeDNA separately, then link them to a Gibson
Assembly station: 

.. raw:: html

    <p align="center">
    <img alt="DNA Weaver Logo" title="DNA Weaver Logo"
         src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaWeaver/master/docs/_static/images/two_vendors_supply_network.png" width="250"/>
    <br /><br />
    </p>

.. code:: python

    import dnaweaver as dw

    cheap_dna_offer = dw.CommercialDnaOffer(
        name="CheapDNA",
        sequence_constraints=[
            dw.NoPatternConstraint(enzyme="BsaI"),
            dw.SequenceLengthConstraint(max_length=4000)
        ],
        pricing=dw.PerBasepairPricing(0.10),
        lead_time=40
    )
    deluxe_dna_offer = dw.CommercialDnaOffer(
        name="DeluxeDNA",
        sequence_constraints=[dw.SequenceLengthConstraint(max_length=3000)],
        pricing=dw.PerBasepairPricing(0.20),
        lead_time=20
    )
    assembly_station = dw.DnaAssemblyStation(
        name="Gibson Assembly Station",
        assembly_method=dw.GibsonAssemblyMethod(
            overhang_selector=dw.TmOverhangSelector(min_tm=55, max_tm=70),
            min_segment_length=500,
            max_segment_length=4000,
            duration=5
        ),
        dna_source=[cheap_dna_offer, deluxe_dna_offer],
        coarse_grain=20
    )
    sequence = dw.random_dna_sequence(10000, seed=123)
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)

    print (quote.assembly_step_summary())

This code prints out an assembly summary showing the source of the
different sequence segments (start, end):

.. code:: bash

    Ordering plan:
      0-1719: From CheapDNA - price 172.80 - lead_time 40.0
      1719-4429: From CheapDNA - price 273.00 - lead_time 40.0
      4429-5318: From DeluxeDNA - price 182.00 - lead_time 20.0
      5318-7359: From CheapDNA - price 206.00 - lead_time 40.0
      7359-10000: From CheapDNA - price 265.00 - lead_time 40.0
    Price: 1098.80, total lead_time: 45.0

Notice how DNA Weaver uses preferentially CheapDNA, with the exception of a 1kb
fragment in the middle of the sequence, which had to be ordered from DeluxeDNA
due to the presence of a BsaI site.

Multi-step assembly
~~~~~~~~~~~~~~~~~~~~~

By defining more DNA sources and connecting them together it is possible to
model complex assembly problems.

For instance in `this example <>`_ we implement a complex DNA assembly chain,
where the final DNA sequence (typically 50kb) is obtained from Yeast
recombination of DNA chunks originating either from the E. coli chromosome
(via PCR extraction) or from the assembly of smaller fragments
via Golden Assembly or Gibson assembly (whichever method is best adapted). These
assembly fragments are obtained either from commercial providers (CheapDNA and
DeluxeDNA) or assembled from oligos:

.. raw:: html

    <p align="center">
    <img alt="DNA Weaver Logo" title="DNA Weaver Logo"
         src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaWeaver/master/docs/_static/images/multiple_steps_supply_network.png" width="250"/>
    <br /><br />
    </p>

Just a few lines of code can produce a comprehensive report (see a sample `here <>`_)
featuring plots of the final assembly plan , comprehensive PDF reports
listing all operations needed, and genbank/fasta files of the sequences to order:

.. code:: python

    from dnaweaver.reports import JsonQuote, make_folder_report
    ...
    quote = assembly_station.get_quote(sequence, with_assembly_plan=True)
    quote.compute_full_assembly_tree()
    json_quote = JsonQuote.from_dnaweaver_quote(quote)
    make_folder_report(json_quote, "report.zip")

Result:

.. raw:: html

    <p align="center">
    <img alt="DNA Weaver Logo" title="DNA Weaver Logo"
         src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaWeaver/master/docs/_static/images/report_illustration.png" width="250"/>
    <br /><br />
    </p>

Assembly with more or less parts reuse
--------------------------------------




Installation
-------------

You can install DnaWeaver through PIP
::
    sudo pip install dnaweaver

Alternatively, you can unzip the sources in a folder and type
::
    sudo python setup.py install

Also install the ncbi-blast+ package to be able to use PCR stations. On Ubuntu:

::
    sudo apt-get install ncbi-blast+

Reports generation needs more dependencies for plots and tables. Install Python dependencies with:

::
    sudo pip install pandas dna_features_viewer weasyprint

You may also need the following non-python dependencies for report generation,
on Ubuntu:

::
    sudo apt-get installbuild-essential python3-dev python3-pip \
        python3-cffi libcairo2 libpango-1.0-0 libpangocairo-1.0-0 \
        libgdk-pixbuf2.0-0 libffi-dev shared-mime-info

License = MIT
--------------

DNA Weaver is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Weaver is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.

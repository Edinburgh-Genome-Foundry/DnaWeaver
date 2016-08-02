.. raw:: html

    <a href="https://twitter.com/share" class="twitter-share-button"
    data-text="DnaWeaver - A Python module for printing with living matter" data-size="large" data-hashtags="Bioprinting">Tweet
    </a>
    <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
    if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
    fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
    </script>
    <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=dnaweaver&type=watch&count=true&size=large"
    allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>

DNA Weaver
============

DnaWeaver is a Python library to determine which DNA fragments to orders from
synthesis company in order to assemble large DNA sequences using Gibson Assembly,
Golden Gate assembly, etc..

Provided a DNA sequence and a list of offers (i.e. pricings and constraints) from
different DNA companies, DnaWeaver returns the list of fragments orders which
minimizes the total cost.

DnaWeaver uses decomposition optimization techniques based on graphs, see
:ref:`howitworks` for more details on how it works.

DnaWeaver can be easily extended to include new offers and constraints from DNA company,
or to take into account objectives other than just price.

Minimal example
---------------
::

    from dnaweaver import *
    from dnachisel import PROVIDERS_CONSTRAINTS
    offer_1 = DnaOrderingProblem("offer 1",
                                 constraints=PROVIDERS_CONSTRAINTS["IDT"],
                                 pricing=lambda sequence: 0.1 * len(sequence))
    offer_2 = DnaOrderingProblem("offer 2",
                                 constraints=PROVIDERS_CONSTRAINTS["gen9"],
                                 pricing=lambda sequence: 0.2 * len(sequence))
    sequence = open("sequence.txt","r").read()
    problem = DnaOrderingProblem(sequence, offers=[offer_1, offer_2],
                                 assembly_method=GibsonAssemblyMethod(20))
    solution = problem.solve(min_segment_length=100, max_segment_length=4000)
    print (solution.summary())

Installation
-------------

You can install DnaWeaver through PIP
::
  sudo pip install DnaWeaver

Alternatively, you can unzip the sources in a folder and type
::
  sudo python setup.py install


Contribute
----------

DnaWeaver is an open-source library originally written at the Edinburgh Genome Foundry by Zulko_.
It is released on Github under the MIT licence, everyone is welcome to contribute.



.. toctree::
    :hidden:
    :maxdepth: 3

    self

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 3

    how_it_works
    ref

.. toctree::
    :caption: Examples
    :maxdepth: 1

    examples/basic_example
    examples/realistic_companies_example
    examples/oligos_assembly_example
    examples/golden_gate_simple_example
    examples/golden_gate_with_forced_cuts_example


.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/DnaWeaver
.. _PYPI: https://pypi.python.org/pypi/DnaWeaver

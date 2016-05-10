.. raw:: html

    <a href="https://twitter.com/share" class="twitter-share-button"
    data-text="DnaAdvisor - A Python module for printing with living matter" data-size="large" data-hashtags="Bioprinting">Tweet
    </a>
    <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
    if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
    fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
    </script>
    <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=dnaadvisor&type=watch&count=true&size=large"
    allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>

DNA Advisor
============

DnaAdvisor is a Python library to determine which DNA fragments to orders from
synthesis company in order to assemble large DNA sequences using Gibson Assembly,
Golden Gate assembly, etc..

Provided a DNA sequence and a list of offers (i.e. pricings and constraints) from
different DNA companies, DnaAdvisor returns the list of fragments orders which
minimizes the total cost.

DnaAdvisor can be easily extended to include new offers and constraints from DNA company,
or to take into account objectives other than just price.

Example of use
---------------



How it works
---------------


Graph representation of the problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaAdvisor uses the following strategies to avoid exploring the whole search space:

.. figure:: images/base_problem.png
   :figwidth: 75%
   :align: center

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaAdvisor uses the following strategies to avoid exploring the whole search space:

.. figure:: images/graph.png
   :figwidth: 75%
   :align: center

Cuts refinement
~~~~~~~~~~~~~~~

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaAdvisor uses the following strategies to avoid exploring the whole search space:
using again the graph trick. The graph of a cuts refinement problem looks like this:

.. figure:: images/refinement.png
   :figwidth: 60%
   :align: center

Installation
-------------

You can install DnaAdvisor through PIP
::
  sudo pip install DnaAdvisor

Alternatively, you can unzip the sources in a folder and type
::
  sudo python setup.py install


Contribute
----------

DnaAdvisor is an open-source library originally written at the Edinburgh Genome Foundry by Zulko_.
It is released on Github under the MIT licence, everyone is welcome to contribute.



.. toctree::
    :hidden:
    :maxdepth: 3

    self

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 3

    ref

.. toctree::
    :caption: Examples

    examples/basic_example
    examples/real_companies_example
    examples/melting_temperature


.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/DnaAdvisor
.. _PYPI: https://pypi.python.org/pypi/DnaAdvisor

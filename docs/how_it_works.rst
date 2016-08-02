.. _howitworks:

How DNA Weaver works
----------------------

Solving a DNA ordering problem using graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaWeaver uses the following strategies to avoid exploring the whole search space:

.. figure:: images/base_problem.png
   :figwidth: 75%
   :align: center

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaWeaver uses the following strategies to avoid exploring the whole search space:

.. figure:: images/graph.png
   :figwidth: 75%
   :align: center

DNA ordering problem for nested assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(All figures to check later with Chantal.)

- The final sequence (60kb) will be done by assembling ~10

.. figure:: images/hierarchical_assemblies.png
   :figwidth: 95%
   :align: center

Each assembly station solves its own DNA ordering problems using its own graphs.
The weight of each edges of the last assembly station's graph is computed by
finding a shortest path in the graph of the second-to-last assembly station.
The second-to-last station gets its edge weights by solving the third-to-last graph problem, etc.

.. figure:: images/hierarchical_assemblies_graphs.png
  :figwidth: 80%
  :align: center

While running this method, the different assembly stations may be presented many
times with the same sequences to assemble optimally. Therefore we use hash-tables to
quickly retrieve the price of already-evaluated sequences.


Cuts refinement
~~~~~~~~~~~~~~~

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
Dna Weaver uses the following strategies to avoid exploring the whole search space:
using again the graph trick. The graph of a cuts refinement problem looks like this:

.. figure:: images/refinement.png
  :figwidth: 60%
  :align: center

Penalty on the number of segments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaWeaver uses the following strategies to avoid exploring the whole search space:
using again the graph trick. The graph of a cuts refinement problem looks like this:

.. figure:: images/refinement.png
  :figwidth: 60%
  :align: center


 (figure with the effect of increassing the penalty)

Forced cuts
~~~~~~~~~~~~~
It is possible to force cuts at certain positions.

(figure)

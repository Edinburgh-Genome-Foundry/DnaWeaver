DNA Weaver assembly guide
--------------------------

Thanks for using DNA Weaver ! This directory contains everything to guide you
through your assembly.

Below is the list of files featured in the directory:


Assembly_report.pdf
```````````````````
This report gives general informations on the sequence and the assembly, as well
as the detail of every assembly step.


sequences.csv
`````````````
This spreadsheet gives the sequence of every fragment and intermediary construct
in the assembly. It can be used to order all basic fragments by copy/pasting the
sequences into an ordering spreadsheet.


genbank/
````````
Contains every fragment and intermediary construct in a genbank format to be
easily browsable using a sequence viewer, and allow for alignment with
sequencing results.


figures/assembly_graph.pdf
``````````````````````````
Complete tree of the assembly, featuring all assembly steps and their
dependencies.


figures/assembly_blocs.pdf
``````````````````````````
Similar to the assembly graphs, but the intermediary constructs are represented
using blocs which give an indication of their respective sizes.


figures/assembly_timeline.pdf
``````````````````````````
If a deadline was supplied for the assembly, prints a Gantt chart of the
assembly, showing the deadline for each intermediary construct. (still an
experimental feature at this stage)


figures/supply_network.pdf
``````````````````````````
Schematic representation of the supply network of the problem

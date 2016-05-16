.. realistic_companies_example:

A example with realistic companies constraints
-----------------------------------------------

In this example we show how to use DnaChisel to model the constraints of companies
when defining DNA offers.

We want to order fragments to assemble a 700bp piece of DNA, and we consider two
fictitious companies which have many complex rules like in real life:

**Company Cell9** charges 0.10$/bp and applies the following constraints:

- GC content must be between 40-65% (and 25-80%, over 50bp windows)
- No BsaI or AarI sites.
- No homopolymers of size 8+ (for A,T,C) or 6+ (for G).

**Company TDI** charges 0.20$/bp and applies the following constraints:

- GC content must be 40-68% globally, 28-76% over 100bp windows, 15-90%
  over 20bp windows, and 24-76% in the 30bp terminal windows.
- No 3-mers repeated 5+ times, or 2-mers repeated 9+ times
- No 6+ homopolymers of G-C or 9+ homopolymers of A-T
- No hairpins (defined as 20-mers with a reverse-complement in the next 200bp window)
- No homopolymers of size 8+ (for A,T,C) or 6+ (for G).

Code
~~~~

.. literalinclude:: ../../examples/realistic_companies_example.py

Expected output:

(Fast) graph homomorphism and core searching
============================================

Works with graphs from Sage and NetworkX.
Also works with plain adjacency lists of graphs on vertices 0,1,... n-1 (C++ interface).

The core branching alg. is implemented in C++, with Cython interface to Python.
Accepts partial maps. Can find both homomorphisms and retracts (with special speedups).


Heuristics include second neighborhood candidates Optionally looks for homomorphisms in the target graph.

Includes heuristics for branching vertex ordering, pruning map candidates based on first and second neighborhood
of a mapped vertex, more heuristics for retract search.

Usage
-----

Run `make` to compile the python module.

Run `make test` to compile and run both the C++ and Python testsuites.

Import with `import homsearch`, use `homsearch.find_homomorphisms` and `homsearch.find_retracts`.

License
-------

This work is licensed under [The MIT licese](http://opensource.org/licenses/MIT).

If you find this code to be of any use, please let me know. I would also welcome any feedback.

Tomas Gavenciak, gavento@ucw.cz

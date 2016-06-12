# flash
Simple C++ implementation of a numerical lattice theory [1,2] for polymer adsorption taking into account the weight distribution of (co)polymers.

![example](/pics/example.png)

Binary for Linux 64bit provided. Run with:

./flash

Modify files SETUP (interaction parameters) and distribution (polymer sequence)

Plot results with gnuplot. E.g. 

gnuplot ps_ens3d.plt

To compile the code use ./build

[1] Statistical theory of the adsorption of interacting chain molecules. 1. Partition function, segment density distribution, and adsorption isotherms
J. M. H. M. Scheutjens and G. J. Fleer
The Journal of Physical Chemistry 1979 83 (12), 1619-1635

[2] Polymers at interfaces, by G. J. Fleer, M. A. Cohen Stuart, J. M. H. M. Scheutjens, T. Cosgrove and B. Vincent. Chapman and Hall, London, 1993. Pp. xv + 502. ISBN 0-412-58160-4




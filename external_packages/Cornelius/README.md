# Cornelius

This version of the Cornelius algorithm is an improved and tested version 
taken from this repository: [Hendrik1704/Cornelius](https://github.com/Hendrik1704/Cornelius)
Commit hash: `93e2c2cf7db39d1578f3794729e24818e2478bf9`

This code is based on cornelius++ version 1.3 (Copyright 2012) by Pasi Huovinen and Hannu Holopainen.
The Cornelius algorithm is aimed to be used as a part of the fluid dynamical models of the heavy-ion physics community.
Publications using this algorithm should cite:
 * P. Huovinen and H. Petersen, Eur.Phys.J.A 48 (2012) 171, arXiv:1206.3371 [nucl-th]

This new version leverages new C++17 features and introduces a few unit tests.
The main change is the one from manual dynamic memory management to automatic 
dynamic memory management in some parts of the code using `std::vector`, while 
most arrays are implemented as `std::array` now.
There are also tests in the separate GitHub repository making sure that the results obtained by the new version 
are compatible with the old one.

Modifications to the original code are done by:
 * Hendrik Roch
 * Haydar Mehryar
 * Joe Latessa
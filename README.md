# CRNcode

The code in this package was developed between 2010 and 2022 and reflects many eras and projects. Most is not carefully error checked. The main focus is the theory of chemical reaction networks (CRNs) - although there is some code whose utility goes beyond CRNs, most was originally inspired by this application. 

###Prerequisites

It relies on a number of C/C++/linux packages and utilities. These include 

1) [GiNaC for symbolic algebra](https://ginac.de/) which I installed on Ubuntu using:
`$ sudo apt install libginac-dev*`

2) [GLPK for linear programming](https://www.gnu.org/software/glpk/) which I installed on Ubuntu using:
`$ sudo apt install glpk-utils libglpk-dev glpk-doc python-glpk`

3) [csdp for semidefinite programming](https://github.com/coin-or/Csdp) which I installed on Ubuntu using:
`$ sudo apt install coinor-csdp`

4) [NAUTY for graph isomorphism](https://users.cecs.anu.edu.au/~bdm/nauty/) which I installed on Ubuntu using:
`$ sudo apt install nauty`

csdp and NAUTY and not integrated into the code: they are called using system calls. (In fact there is an over-reliance on system calls in many parts of the code.)

###Creating the binaries

If all goes well, the main binaries can be constructed with:

`$ make`

###Documentation

The [manual](https://github.com/CRNcode/CRN/blob/main/docs/CRNcode.pdf) describes some features of the library. It is currently incomplete. 

###Folders associated with particular projects

The [HopfAtoms folder](https://github.com/CRNcode/CRN/tree/main/HopfAtoms) contains scripts and data associated with the preprint *The smallest bimolecular mass action reaction networks admitting Andronov-Hopf bifurcation* (<https://arxiv.org/abs/2207.04971>) (with Balázs Boros; other scripts associated with this paper are on GitHub [here](https://github.com/balazsboros/reaction_networks/tree/main/3species_4reactions)).


# FHEZ: Bootstrapping FHE over Z (the integers) in less than one second

Proof-of-concept implementation of the fast bootstrapping for fully homomorphic encryption (FHE) schemes over the integers described in the paper [Per21].

The branch PKC21 corresponds to the code used to collect the running times and memory usage shown in this paper, thus, to garantee the reproducibility of the experiments, this branch should not be changed and all possible future modifications will be added to the master branch.

The main classes of this project are:
	* ScalarNandHE - Base scheme, the one that is bootstrapped. It can perform one NAND gate before bootstrapping.
	* GAHE - Scheme used to bootstrap the base scheme.
	* BootstrapperSingleNandHE - Implements the bootstrapping.

Dependencies: install NTL (https://libntl.org/) and GMP with C++ support enabled (search for --enable-cxx in https://gmplib.org/manual/Build-Options).

To compile, just run make, which will create three executable files with self-explanatory names.

## References:

[Per21] Hilder Vitor Lima Pereira. Bootstrapping fully homomorphic encryption over the integers in less than one second. Published in [PKC 2021](https://pkc.iacr.org/2021/).

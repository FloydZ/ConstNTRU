Goal
-----
The goal of this repo is to implement a `C++ constexpr` version of the NIST round 3 submission [NTRU](https://ntru.org/).

WARNING
-----
THIS IS NOT SUITED FOR REAL APPLICATIONS. DO NOT USE IT IN PRODUCTION. You have been warned. 

Setup
----
clone the repository via:
```bash
 git clone --recursive https://github.com/FloydZ/ConstNTRU
```
Note that currently the implementation depents on [ConstRand] (https://github.com/FloydZ/ConstRand) repository, which 
exports the function 
```
auto rand = const_uniform_distribution_array<int32_t, 2*N>(-(Q/2), Q/2)
```
This function is important, because the NTRU implementation needs an array of random byts which is are generated at compile time.

If you dont want you do not want to depent on this repository, you can replace
this funciton by any other wich outputs an array of type `int32_t[2*N]`.

Parameters
----
The following 4 parameter sets are implemented
    - NTRUHPS2048509
    - NTRUHPS2048677
    - NTRUHPS4096821
    - NTRUHRSS701
You can select on of them via the preprocessor flag `-D ...` or `#define ...`.

FPLLL
---
If you define the preprocessor flag `FPLLL_H` (which is defined in `fplll.h` in the [fplll](https://github.com/fplll/fplll) repository)
additonally functions 
```bash
// create a NTRU matrix from a pk `h`
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &m, const int32_t h[N]);

// create a new NTRU matrix from the seed 0
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &A);

// create a new NTRU matrix from the seed 0
void generate_const_ntru_encrypt_matrix_const(fplll::ZZ_mat<ZT> &A);
```
getting available. 

TODO
----
- implement KAT test from [here](https://ntru.org/index.shtml)
- explain whats `N` and `Q` are.

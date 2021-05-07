Goal
-----
The goal of this repo is to implement a `C++ constexpr` version of the NIST round 3 submission [NTRU](https://ntru.org/).

WARNING
-----
THIS IS NOT SUITED FOR REAL APPLICATIONS. DO NOT USE IT IN PRODUCTION. You have been warned. 
Maybe you receive warning regarding an integer overflow in `deps/ConstRand/src/const_rand.h`. Ignore it, thts just the rng.

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
This is important, because this NTRU implementation needs an array of random bits which is are generated at compile time.
If you dont like this code, you can replace it with whatever function you would like, the only important point, is that
pass a const `int32_t[2*N]` array.

Parameters
----
The following 4 parameter sets are implemented
    - NTRUHPS2048509
    - NTRUHPS2048677
    - NTRUHPS4096821
    - NTRUHRSS701
You can set the preprocessor flag via `-D ...` or `#define ...`.

FPLLL
---
If you define the preprocessor flag `FPLLL_H` (which is defined in `fplll.h` of the [fplll](https://github.com/fplll/fplll) repository)
additonally functions 
```bash
// create a NTRU matrix from a pk `h`
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &m, const int32_t h[N]);

// create a new NTRU matrix from the seed 0
void generate_const_ntru_encrypt_matrix(fplll::ZZ_mat<ZT> &A);

// create a new NTRU matrix from the seed 0
void generate_const_ntru_encrypt_matrix_const(fplll::ZZ_mat<ZT> &A);
```

TODO
----
- implement KAT test from [here](https://ntru.org/index.shtml)
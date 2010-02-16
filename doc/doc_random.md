% RNP Documentation -- Random
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# RNP::Random

The Random module contains random number generators.
All functions are in namespace `RNP::IO`.

## Contents

* [RandomRealsUniform01](#Random_RandomRealsUniform01)

---
## RandomRealsUniform01<a name="Random_RandomRealsUniform01" />

Generates a random number in the open interval (0,1);

This function corresponds to LAPACK function `dlaruv`.
Here are the LAPACK notes:

This routine uses a multiplicative congruential method with modulus
2^48 and multiplier 33952834046453 (see G.S.Fishman,
'Multiplicative congruential random number generators with modulus
2^b : an exhaustive analysis for b = 32 and a partial analysis for
b = 48', Math. Comp. 189, pp 331-344, 1990).

48-bit integers are stored in 4 integer array elements with 12 bits
per element. Hence the routine is portable across machines with
integers of 32 bits or more.

### Prototype

	template <class T>
	void RandomRealsUniform01(size_t n, T *x, int iseed[4] = NULL);

### Arguments:

=n=
	The number of random numbers to be generated.
=x=
	The array where the generated numbers are to be stored.
=iseed=
	Each entry of the seed must be between 0 and 4095 (inclusive), and `iseed[3]` must be odd.
	On exit, the seed is updated.
	If no seed is passed in (passing in `NULL`), then an internal seed is used and updated.


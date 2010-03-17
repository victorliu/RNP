% RNP -- Rapid Numerical Prototyping
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# RNP -- Rapid Numerical Prototyping

A simple library for common numerical computations.

The [design rationale](about.html) behind the library was to make it maximally easy to use, with performance a secondary consideration.
It turns out however, that for most high level algorithms, there is only one good way to implement them, so the performance tends to be quite good.

## Status

Currently, RNP supports all common linear algebra operations typically used in the physical sciences.
In this respect, RNP is currently feature complete (the target feature set does not include the TODO items marked below), but some interfaces are not finalized, certain parameters are not respected, and memory allocation is not optimal.

## Documentation

Documentation for the different modules is listed below.

* [IO](doc_io.html)
* [Random](doc_random.html)
* [TBLAS](doc_tblas.html)
* [Sparse](doc_sparse.html)
* [LinearSolve](doc_linsolve.html)
* [Eigensystems](doc_zge.html)
* [Generalized Eigensystems](doc_zgg.html)
* [Iterative Linear Solvers](doc_itersolve.html)
* [Iterative Generalized Eigensystems](doc_iterzgg.html)
* [LeastSquares](doc_leastsquares.html) - TODO: Dense least squares/norm approximation
* [IRA](doc_ira.html) (Implicitly Restarted Arnoldi; ARPACK)
* SVD - TODO: Dense SVD
* Iterative SVD - TODO
* Integration - TODO

Documentation for activities related to the library is listed below.

* [Compilation and file dependencies](doc_compiling.html)
* [Porting to heavy duty numerical packages](doc_porting.html)

## Download

* From now on, the code is [hosted on GitHub](http://github.com/victorliu/RNP). You can download the latest snapshot there.
* [rnp-20100116](rnp_20100116.zip) - Initial release.

## FAQ

* Are there specialized routines for symmetric or positive definite problems?
  
  No, not yet. Chances are you can find code on the Internet for symmetric routines, since those are far simpler.

* How do I invert a matrix?
  
  Perform a matrix solve with your matrix and use the identity matrix for the right hand side.
  A slightly more optimal way is to perform the factorization and invert the L and U factors yourself.
  Seriously, this is the best way. There is no built function since matrix solves are destructive, and the storage conventions depend on the application. You can write your own suitable wrapper for this if you really want it to be more convenient.

## TODO

* Ensure IterativeGeneralizedEigensolvers respect iteration limits and multiply limits.
* Improve tests.
* Port some symmetric eigensolvers: zheev, dsyev, and PRIMME's JDQMR.
* Port root-finding routines.
* Rework JDQZ into its own namespace.
* Add constrained JDQZ (linear homogeneous constraints).

## License

I'm not big on keeping track of licenses, but I will do my best below.
Since this library is intended for prototyping, I think it is unlikely anyone will have to bother with looking at licenses.

The code contained in RNP is a mixture of code I wrote and code translated from other sources.
Modules IO, Sparse, and the Jacobi eigenssolver are written by myself, and are in the public domain.
Modules Random, most of TBLAS, LinearSolve, Eigensystems, and GeneralizedEigensystems are translated from BLAS/LAPACK, and are subject to its license. IterativeLinearSolvers and IterativeGeneralizedEigensystems are translated from Matlab/Fortran code by Martin van Gijzen, Peter Sonneveld, and Gerard L. G. Sleijpen, and are under the GNU General Public License version 2 or later.

## Author/collector

This library is maintained by Victor Liu. Email is vee kay ell (three letters) at stanford dot edu.

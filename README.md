# RNP -- Rapid Numerical Prototyping

A simple library for common numerical computations.

The [design rationale](about.html) behind the library was to make it maximally easy to use, with performance a secondary consideration.
It turns out however, that for most high level algorithms, there is only one good way to implement them, so the performance tends to be quite good.

## Status

Currently, RNP supports all common linear algebra operations typically used in the physical sciences.
In this respect, RNP is currently feature complete (the target feature set does not include the TODO items marked below), but some interfaces are not finalized, certain parameters are not respected, and memory allocation is not optimal.

## Documentation

Documentation for the different modules is listed below. (Note that GitHub does not support the Discount Markdown extensions, so some of these pages will not render properly; view the actual HTML files for actual reference).

* [IO](http://github.com/victorliu/RNP/tree/master/doc/doc_io.md)
* [Random](http://github.com/victorliu/RNP/tree/master/doc/doc_random.md)
* [TBLAS](http://github.com/victorliu/RNP/tree/master/doc/doc_tblas.md)
* [Sparse](http://github.com/victorliu/RNP/tree/master/doc/doc_sparse.md)
* [LinearSolve](http://github.com/victorliu/RNP/tree/master/doc/doc_linsolve.md)
* [Eigensystems](http://github.com/victorliu/RNP/tree/master/doc/doc_zge.md)
* [Generalized Eigensystems](http://github.com/victorliu/RNP/tree/master/doc/doc_zgg.md)
* [Iterative Linear Solvers](http://github.com/victorliu/RNP/tree/master/doc/doc_itersolve.md)
* [Iterative Generalized Eigensystems](http://github.com/victorliu/RNP/tree/master/doc/doc_iterzgg.md)
* [LeastSquares](http://github.com/victorliu/RNP/tree/master/doc/doc_leastsquares.md) - TODO: Dense least squares/norm approximation
* SVD - TODO: Dense SVD
* Iterative SVD - TODO
* Integration - TODO

Documentation for activities related to the library is listed below.

* [Compilation and file dependencies](http://github.com/victorliu/RNP/tree/master/doc/doc_compiling.md)
* [Porting to heavy duty numerical packages](http://github.com/victorliu/RNP/tree/master/doc/doc_porting.md)

## Download

* From now on, the code is [hosted on GitHub](http://github.com/victorliu/RNP). You can download the latest snapshot there.

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
* Rework JDQZ into its own namespace.
* Add constrained JDQZ (linear homogeneous constraints).
* Port some Hermitian LAPACK routines.

## License

I'm not big on keeping track of licenses, but I will do my best below.
Since this library is intended for prototyping, I think it is unlikely anyone will have to bother with looking at licenses.

The code contained in RNP is a mixture of code I wrote and code translated from other sources.
Modules IO, Sparse, and the Jacobi eigenssolver are written by myself, and are in the public domain.
Modules Random, most of TBLAS, LinearSolve, Eigensystems, and GeneralizedEigensystems are translated from BLAS/LAPACK, and are subject to its license. IterativeLinearSolvers and IterativeGeneralizedEigensystems are translated from Matlab/Fortran code by Martin van Gijzen, Peter Sonneveld, and Gerard L. G. Sleijpen, and are under the GNU General Public License version 2 or later.

## Author/collector

This library is maintained by Victor Liu. Email is vee kay ell (three letters) at stanford dot edu.

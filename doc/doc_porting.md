% RNP Documentation -- Porting to heavy duty numerical packages
% Victor Liu (vkl@stanford.edu)
% February 12, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Porting to heavy duty numerical packages

Many of the routines are taken directly from Fortran packages, with the interfaces kept fairly close to the original, so porting to an industrial strength, optimized numerical library should be straightforward.
Below is listed the original libraries from which each module originated.

In many instances, simply optimizing the TBLAS layer should afford substantial improvements in performance.
This can be done relatively simply by providing explicit template specializations of TBLAS functions.
Doing this will not achieve performance on par with optimized LAPACK libraries since all blocking has been stripped from the LAPACK routines.
Future developments may restore some blocking functionality.

The BLAS and LAPACK bindings to the core API are being added right now.
In the near future, using vendor provided versions of these libraries will be as simple as defining the preprocessor macros `RNP_HAVE_BLAS` and `RNP_HAVE_LAPACK`.
Support currently exists for most of level 1 BLAS, and complex `zgemv` and `zgemm` for non-transposed matrices.

## Library origins

* IO - None.
* Random - From LAPACK. Any respectable pseudorandom number generator should be easily substitutable.
* TBLAS - From official BLAS. Note that all interfaces are identical to the originals, except for character job options, which were converted to template parameters in the same order as the original.
* Sparse - None. The matrix format is generic compressed row or column format, so most sparse matrix libraries can handle this.
* LinearSolve - From LAPACK. Follows the LAPACK interface closely.
* Eigensystems - From LAPACK; a direct translation, but removing the character job switches and info return. Note that Eigensystems_jacobi follows the LAPACK interface, but requires no workspace.
* Generalized Eigensystems - From LAPACK; follows the LAPACK interface like Eigensystems.
* Iterative Linear Solvers - Due to the diversity of implementation and interfaces of iterative linear solvers, porting these are not very straightforward. Most iterative linear solvers support the same set of options, so it would simply be a matter of rewriting the function calls. In particular, most solvers use the convention that the matrix multiplication routine returns the result into a new vector, while the preconditioner does not, which is the convention used in RNP.
* Iterative Generalized Eigensystems - To the author's knowledge, there are only a handful of methods that can deal with the fully general case, and none of them allow for arbitrary eigenvalue sorting like RNP. Those that come close are ARPACK (or ARPACK++ for C++) and JDQZ (which is the method from which the RNP implementation was derived). The Matlab JDQZ supports many more options than what is implemented here, but it is unlikely that anyone prototypes in C++ only to use Matlab in the end. In that respect, they are probably stuck with RNP. The interface to ARPACK++ does not require a linear solver, and only certain predefined eigenvalue targets are permitted. Also, preconditioning is not handled explicitly, so substantial consultation of the documentation may be required.

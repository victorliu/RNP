% RNP -- About
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Rapid Numerical Prototyping

Ever wanted to try out some numerical algorithm, but dread having to link your code with some archaic Fortran libraries?
The Rapid Numerical Prototyping library aims to solve this problem by making the easiest-to-use numerical C++ library.
RNP implements algorithms in a straightforward way when possible, and always chooses the most generic or simplest implementation.
All implementation/source files are self-contained, so a minimal number of files ever need to be included.

## Design rationale

Here I will explain the reasons for the API/interfaces for various components.

### Dense matrix operations and BLAS implementation

* Operator overloading is not used since the C++ operator set is insufficiently rich to express a large enough subset of matrix operations.
  For example, multiplication by (conjugate) transpose cannot be directly expressed using overloadable operators.
  Attempts to get around this shortcoming with templates, matrix expression objects, etc. end up producing hard-to-understand template code that is equally hard to debug.
* All the character arguments (`job` or `uplo`) are template parameters since I have never seen code that needed to choose them at runtime.
* Only operations on general (no symmetry) matrices are supported for simplicity.
* BLAS implementations provide no argument checking, perform no special casing, and are usually about 3-5 lines of code to try to offset some of the performance lost to using unblocked algorithms.
* Matrix formats match Fortran interface expectations so that prototypes built with RNP are easily portable to optimized Fortran libraries.

### Dense eigensolvers

* The dense eigensystem solver (equivalent of `zgeev`) implements a Jacobi-like algorithm, which is amazing simple.
  The standard LAPACK implementation performs QR iteration, but that code contains many oddities (such as manually unrolled recursion).
  The Jacobi-like algorithm will probably produce more accurate eigenvalues/vectors, at the cost of being somewhat slower for large problems.
* The dense generalized eigensystem solver (equivalent of `zggev` or `zgevv`) is a direct translation of the LAPACK QZ iteration implementation.
The code has been converted using a modified F2C and hand re-written to remove all goto's and operate out of a single source file.

### Linear solvers

* The dense linear solver is a direct LU decomposition with pivoting and triangular solve, like in LAPACK.
* Several (sparse) iterative linear solvers are implemented: IDR(s), BiCGSTAB(l), and GMRES. These are the current state-of-the-art non-symmetric iterative linear solvers.
* There are no plans to implement any direct solvers for large sparse systems due to the overwhelmin complexity of those algorithms.

### Sparse eigensolvers

* There is only a single eigensolver implemented: the JDQZ algorithm for generalized non-symmetric eigensystems.
No non-generalized eigensolver is available due to the desire to support arbitrary spectral selection (you can choose eigenvalues using a user-defined eigenvalue sorter); use the generalized algorithm with the B matrix equal to the identity if needed.
The sparse eigensolver takes function pointers to perform matrix-vector multiplies, linear solves, and preconditioning.
All of these operations can be performed with other routines in RNP, and adapters are provided for convenience.

### Dense SVD

* I plan to port the LAPACK Jacobi SVD algorithm for this purpose.


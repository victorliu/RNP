% RNP Documentation -- Eigensystems
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Eigensystems

This module contains functions related to the computation of eigensystems or Schur decompositions.
All functions are in namespace `RNP`.

## Contents

* [Eigensystem](#RNP_Eigensystem)
* [Eigensystem_jacobi](#RNP_Eigensystem_jacobi)

---
## Eigensystem<a name="RNP_Eigensystem" />

Eigensystem computes for an `n`x`n` complex nonsymmetric matrix A, the
eigenvalues and, optionally, the left and/or right eigenvectors.

The right eigenvector `v[j]` of A satisfies
                 `A * v[j] = lambda[j] * v[j]`
where `lambda[j]` is its eigenvalue.
The left eigenvector `u[j]` of A satisfies
              `u[j]^H * A = lambda[j] * u[j]^H`
where `u[j]^H` denotes the conjugate transpose of `u[j]`.

The computed eigenvectors are normalized to have Euclidean norm
equal to 1 and largest component real.

This function corresponds to LAPACK functions `zgeev` and `cgeev`.

Return value:

* `= 0`:  Successful exit
* `< 0`:  if return value = `-i`, the (1-based) `i`-th argument had an illegal value.
* `> 0`:  the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed; elements and `return value`,...,`n-1` of eval contain eigenvalues which have converged.

### Prototype

	int Eigensystem(
		size_t n, std::complex<double>* a, size_t lda,
		std::complex<double>* eval,
		std::complex<double>* evec, size_t ldvec,
		std::complex<double> *work = NULL, double *rwork = NULL,
		size_t lwork = 0);

### Arguments:

=n=
	The number of rows/columsn of matrix `A`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=eval=
	Pointer to output space for eigenvalues, stored contiguously.
=vl=
	Pointer to output matrix for left eigenvectors, stored in columns in same order as `eval`.
	If `NULL`, then no left eigenvectors will be computed.
=ldvl=
	Leading dimension of matrix `vl`.
=vr=
	Pointer to output matrix for right eigenvectors, stored in columns in same order as `eval`.
	If `NULL`, then no right eigenvectors will be computed.
=ldvr=
	Leading dimension of matrix `vr`.
=work=
	Pointer to complex workspace, should be length `2*n`.
	If `NULL`, the workspace is internally allocated and freed.
=rwork=
	Pointer to real workspace, should be length `2*n`.
	If `NULL`, the workspace is internally allocated and freed.
=lwork=
	Length of the workspace parameter `work`. If set to `0`, then it is assumed that `work` is of sufficiently large size.
	If set to `-1`, a workspace query is performed, and the function returns immediately and the optimal workspace size for the provided arguments is stored in `work[0]`.
	The other pointer arguments are not referenced.

## Eigensystem_jacobi<a name="RNP_Eigensystem_jacobi" />

This function is identical to [Eigensystem](#RNP_Eigensystem) except that it uses modified Jacobi rotations to diagonalize the matrix and requires no workspace.
This is substantially slower than the QR iteration of the standard LAPACK algorithm above, but it may provide better relative accuracy of results.

Return value:

* `= 0`: Successful exit.
* `> 0`: Failure to fully diagonalize matrix. Return value is the number of Jacobi passes performed so far.

### Prototype

	int Eigensystem_jacobi(
		size_t n, std::complex<double>* a, size_t lda,
		std::complex<double>* eval,
		std::complex<double>* evec, size_t ldvec,
		std::complex<double> *work, double *rwork);

### Arguments:

=n=
	The number of rows/columsn of matrix `A`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=eval=
	Pointer to output space for eigenvalues, stored contiguously.
=vl=
	Pointer to output matrix for left eigenvectors, stored in columns in same order as `eval`.
	If `NULL`, then no left eigenvectors will be computed.
=ldvl=
	Leading dimension of matrix `vl`.
=vr=
	Pointer to output matrix for right eigenvectors, stored in columns in same order as `eval`.
	If `NULL`, then no right eigenvectors will be computed.
=ldvr=
	Leading dimension of matrix `vr`.
=work=
	This argument is not used.
=rwork=
	This argument is not used.

% RNP Documentation -- Generalized eigensystems
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Generalized eigensystems

This module contains functions related to the computation of generalized eigensystems or generalized Schur decompositions.
All functions are in namespace `RNP`.

## Contents

* [GeneralizedEigensystem](#RNP_GeneralizedEigensystem)
* [GeneralizedSchurDecomposition](#RNP_GeneralizedSchurDecomposition)

---
## GeneralizedEigensystem<a name="RNP_GeneralizedEigensystem" />

Computes for a pair of `n`x`n` matrices (A,B), the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.

A generalized eigenvalue for a pair of matrices (A,B) is a scalar
`lambda` or a ratio `alpha/beta = lambda`, such that `A - lambda * B` is
singular. It is usually represented as the pair (`alpha`,`beta`), as
there is a reasonable interpretation for `beta = 0`, and even for both
being zero.

The right generalized eigenvector `v[j]` corresponding to the
generalized eigenvalue `lambda[j]` of (A,B) satisfies
             `A * v[j] = lambda[j] * B * v[j]`.
The left generalized eigenvector `u[j]` corresponding to the
generalized eigenvalues `lambda[j]` of (A,B) satisfies
             `u[j]^H * A = lambda[j] * u[j]^H * B`
where `u[j]^H` is the conjugate-transpose of `u[j]`.
We denote by `U` the matrix whose columns are `u[j]` for `j = 1 ... N` (and correspondingly for `V`).

This function corresponds to LAPACK functions `zggev`, `zgegv`, `cggev`, and `cgegv`.

Return value:

* `= 0`: Successful exit
* `< 0`: if return value = `-i`, the (1-based) `i`-th argument had an illegal value.
* `=1,...,n+1`: The QZ iteration failed.  No eigenvectors have been calculated, but `alpha[j]` and `beta[j]` should be correct for `j = `return value-1`,...,n-1`.
* `=n+2`: There was an error computing the eigenvectors.

### Prototype

	int GeneralizedEigensystem(size_t n, 
		std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, 
		std::complex<double> *alpha, std::complex<double> *beta,
		std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
		std::complex<double> *work, double *rwork);

### Arguments:

=n=
	The number of rows/columns of matrices `A`, `B`, `The order of the matrices A, B, `U`, and `V`.
=a=
	Pointer to the first element of matrix `A`.
	On exit, a is overwritten.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=b=
	Pointer to the first element of matrix `B`.
	On exit, b is overwritten.
=ldb=
	Leading dimension of matrix `B`.
=alpha, beta=
	Vectors of length `n`.
	On exit, `alpha[j]/beta[j]` will be the generalized eigenvalues.

	Note: the quotients `alpha[j]/beta[j]` may easily over- or
	underflow, and `beta[j]` may even be zero.  Thus, the user
	should avoid naively computing the ratio `alpha/beta`.
	However, alpha will be always less than and usually
	comparable with `norm(A)` in magnitude, and `beta` always less
	than and usually comparable with `norm(B)`.

=vl=
	Pointer to the first element of matrix `U`.
	If `vl != NULL`, the left generalized eigenvectors `u[j]` are
	stored one after another in the columns of `vl`, in the same order as their eigenvalues.
	Each eigenvector is scaled so the largest component has	`abs(real part) + abs(imag part) = 1`.
=ldvl=
	Leading dimension of matrix `U`.
=vr=
	Pointer to the first element of matrix `V`.
	If `vr != NULL`, the right generalized eigenvectors `v[j]` are
	stored one after another in the columns of `vr`, in the same
	order as their eigenvalues.
	Each eigenvector is scaled so the largest component has
	`abs(real part) + abs(imag. part) = 1`.
=ldvr=
	Leading dimension of matrix `V`.
=work=
	Pointer to complex workspace, should be length `2*n`.
=rwork=
	Pointer to real workspace, should be length `8*n`.


---
## GeneralizedSchurDecomposition<a name="RNP_GeneralizedSchurDecomposition" />

Computes for a pair of `n`x`n` matrices (A,B),
the generalized eigenvalues, the generalized complex Schur form
(S, T), and optionally left and/or right Schur vectors (`vsl` and
`vsr`). This gives the generalized Schur factorization
        `(A,B) = ( U*S*V^H, U*T*V^H )`
where we donote by `U` the square matrix whose columns are the left Schur vectors, `V` the square matrix whose columns are the right Schur vectors, and `V^H` is the conjugate-transpose of `V`.

A generalized eigenvalue for a pair of matrices (A,B) is a scalar `w`
or a ratio `alpha/beta = w`, such that  `A - w*B` is singular.  It is
usually represented as the pair (`alpha`,`beta`), as there is a
reasonable interpretation for `beta = 0`, and even for both being zero.

A pair of matrices (S,T) is in generalized complex Schur form if `S`
and `T` are upper triangular and, in addition, the diagonal elements
of `T` are non-negative real numbers.

This function corresponds to LAPACK functions `zgges`, `zgegs`, `cgges`, and `cgegs`.

Return value:

* `= 0`: Successful exit
* `< 0`: if return value = `-i`, the (1-based) `i`-th argument had an illegal value.
* `=1,...,n+1`: The QZ iteration failed.  No eigenvectors have been calculated, but `alpha[j]` and `beta[j]` should be correct for `j = `return value`,...,n-1`.
* `=n+1`: The QZ iteration failed.

### Prototype

	int GeneralizedSchurDecomposition(size_t n,
		std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb,
		std::complex<double> *alpha, std::complex<double> *beta,
		std::complex<double> *vsl, size_t ldvsl, std::complex<double> *vsr, size_t ldvsr,
		std::complex<double> *work, double *rwork);

### Arguments:

=n=
	The number of rows/columns of matrices `A`, `B`, `The order of the matrices A, B, `U`, and `V`.
=a=
	Pointer to the first element of matrix `A`.
	On exit, a is overwritten with matrix `S`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=b=
	Pointer to the first element of matrix `B`.
	On exit, b is overwritten with matrix `T`.
=ldb=
	Leading dimension of matrix `B`.
=alpha, beta=
	Vectors of length `n`.
	On exit, `alpha[j]/beta[j]` will be the diagonals of the complex Schur form (A,B).
	The `beta[j]` will be non-negative real.

	Note: the quotients `alpha[j]/beta[j]` may easily over- or
	underflow, and `beta[j]` may even be zero.  Thus, the user
	should avoid naively computing the ratio `alpha/beta`.
	However, alpha will be always less than and usually
	comparable with `norm(A)` in magnitude, and `beta` always less
	than and usually comparable with `norm(B)`.

=vl=
	Pointer to the first element of matrix `U`.
	If `vl != NULL`, the left generalized eigenvectors `u[j]` are
	stored one after another in the columns of `vl`, in the same order as their eigenvalues.
	Each eigenvector is scaled so the largest component has	`abs(real part) + abs(imag part) = 1`.
=ldvl=
	Leading dimension of matrix `U`.
=vr=
	Pointer to the first element of matrix `V`.
	If `vr != NULL`, the right generalized eigenvectors `v[j]` are
	stored one after another in the columns of `vr`, in the same
	order as their eigenvalues.
	Each eigenvector is scaled so the largest component has
	`abs(real part) + abs(imag. part) = 1`.
=ldvr=
	Leading dimension of matrix `V`.
=work=
	Pointer to complex workspace, should be length `2*n`.
=rwork=
	Pointer to real workspace, should be length `8*n`.


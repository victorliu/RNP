% RNP Documentation -- LinearSolve
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# LinearSolve

The dense linear solvers are all located in namespace `RNP`.

## Contents

* [LinearSolve](#RNP_LinearSolve)
* LinearSolve - TODO: Performs a linear solve without pivot space
* LinearSolveNice - TODO: Performs a linear solve nondestructively, without pivot space
* Invert - TODO: Inverts a matrix
* InvertNice - TODO: In-place inversion

---
## LinearSolve<a name="RNP_LinearSolve" />

Solves the linear system `op(A) * X = B` for dense, square `op(A)`, where `op` is determined by the template parameter `trans`.
`B` is a matrix of right hand sides.

Internally, an LU decomposition with pivoting is used.
Note that specialized versions of this function may be preferrable for better performance on vectors, or for particular conveniences.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameter corresponds to the original LAPACK parameter `trans`, and must be specified.
Here is an example of how the function should be called:

	RNP::LinearSolve<'N'>(n, 1, a, lda, b, n);

This function corresponds to LAPACK functions `zgesv`, `dgesv`, `cgesv`, and `sgesv`.

### Prototype

	template <char trans>
	struct LinearSolve{
		template <class T>
		LinearSolve(size_t n, size_t nRHS,
		            T *a, size_t lda, T *b, size_t ldb,
		            int *info = NULL, size_t *pivots = NULL);
	};

### Arguments:

=trans=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `trans=N`, `op(A) = A`.
	* If `trans=T`, `op(A) = A'` (transpose of A)
	* If `trans=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=n=
	The number of rows/columns of matrix `A`.
=nRHS=
	The number of right hand side vectors, or number of columns of `B`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=b=
	Pointer to the first element of matrix `B`.
=ldb=
	Leading dimension of matrix `B`.
=info=
	If matrix `A` is singular, the position of the first zero pivot encountered in the LU decomposition is return in `info` if it is non-`NULL`. The position is 1-based.
=pivots=
	If non-`NULL`, must be of length `n`. This is used to store pivoting information.
	If it is ignored (pass in `NULL`), then space is automatically allocated and freed.

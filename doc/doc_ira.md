% RNP Documentation -- IRA
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# IRA (Implicitly Restarted Arnoldi; ARPACK)

This module contains a driver interface for a thread-safe C++ port of ARPACK.
The interface is a bit more subtle than `JDQZ`, but the trustworthiness and reliability of this package is more well established.
All functions are in namespace `RNP::IRA`.

## Contents

* [Params](#IRA_Params)
* [Sort Functions](#IRA_SortFunctions)
* [Workspace](#IRA_Workspace)
* [Regular](#IRA_Regular)
* [ShiftInverse](#IRA_ShiftInverse)

---
## Params<a name="IRA_Params" />

The parameter structure describes general settings for the Arnoldi iteration.

=max_iterations=
	The maximum number of Arnoldi iterations allowed.
	On return from a driver function, the actual number of iterations taken is written here (note that this is currently not the case; this will be implemented in the future).
	The default is 1000.
=tol=
	The tolerance or stopping criteria: the relative accuracy of a Ritz value is considered acceptable if the computed relative error bound for the Ritz value is less than this value.
	The default is 0, which indicates that machine precision should be used.
=no_eigenvectors=
	Set this to true to indicate that eigenvectors should not be computed.
	However, the invariant subspace corresponding to the eigenvalues is still returned, so relatively little savings is gained.
	The default is `false`.

---
## Sort Functions<a name="IRA_SortFunctions" />

The sort order of the eigenvalues determines the part of the spectrum that will be returned.
The IRA module generalizes ARPACK's spectrum selection by allowing the user to pass any comparison function.
However, care must be taken in defining a consistent ordering, and furthermore, the proper selection of a driver routine depending on whether interior or extremal eigenvalues are sought.
Several predefined sort functions corresponding to those provided by ARPACK are listed below.

* LargestMagnitude
* SmallestMagnitude
* LargestRealPart
* SmallestRealPart
* LargestImaginaryPart
* SmallestImaginaryPart

None of the predefined sort functions reference the `data` argument.

If a custom `sorter` function is used, then `sorter(a, b)` should return `false` if `a` is more desirable than `b`.

---
## Workspace<a name="IRA_Workspace" />

The workspace allocation is optional, as described below.
However, for repeated calls to the drivers for problems of identical size, it is more efficient to pre-allocate the workspace.
Functions are provided to allocate the properly sized work arrays.

### Prototype

	bool AllocateWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace);
	void FreeWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace);

The first three arguments correspond to the arguments of the same name passed to the driver routines, and should have the same values.
`AllocateWorkspace` returns true upon success.

---
## Regular<a name="IRA_Regular" />

This is the driver routine for regular mode problems (no spectral transformation).
This is the driver to use for computing extremal eigenvalues (e.g. largest magnitude, smallest real part, etc.).
Both generalized and non-generalized problems are handled here, however the operator callback functions must perform different actions in the two cases, so care must be taken in making sure the action matches the problem being solved.

Solves an eigensystem `A * X = lambda * B * X`, where only matrix-vector multiplication routines are available for `A` and `B` (typically the case for large sparse matrices), and a small number of eigenvalues/vectors are needed.
The case of `B = I` is handled specially.

Returns the number of converged eigenpairs.

### Prototype

	int Regular(
		size_t n, Areg Aop, Bfunc Bop,
		size_t n_wanted, size_t n_arnoldi,
		SortFunc sorter,
		std::complex<double> *w,
		std::complex<double> *v, size_t ldv,
		Params *params = NULL,
		Workspace *workspace = NULL,
		void *Adata = NULL, void *Bdata = NULL, void *SortData = NULL
	);

### Arguments:

=n=
	The number of rows/columns of matrices `A` and `B`.
=Aop, Bop=
	If `Bop` is `NULL`, then the non-generalized problem is solved, and `Aop` must perform
		y <- A*x
	If `Bop` is not `NULL`, then the generalized problem is solved, and `Aop` must perform
		y <- inv(B)*A*x
	and `Bop` must perform `y <- B*x`.
	The third argument is user data, and can be passed in via `Adata` and `Bdata`.
	Note that the prototype for `Aop` differs between driver routines.
=n_wanted=
	The number of desired eigenvalues/vectors to be computed (`n_wanted <= n-2`).
=n_arnoldi=
	The maximum number of Arnoldi vectors to use. It is required that `n_wanted+1 <= n_arnoldi <= n`.
	A reasonable choice if unsure is `2*n_wanted+1`.
=sorter=
	The eigenvalue sort function used for selecting which part of the spectrum to compute.
=w=
	Vector of length `n_wanted+1` to store the eigenvalues.
	Only the first `n_wanted` entries should be trusted.
=v=
	Pointer to the first element of the matrix which will hold the eigenvectors.
	It must have at least `n` rows and `n_arnoldi` columns. The eigenvectors are stored in the columns in the same order as the eigenvalues.
	Only the first `n_wanted` vectors should be trusted.
=ldv=
	Leading dimension of matrix `v` (typically `n` unless `v` is a submatrix of a larger matrix).
=params=
	Specifies the IRA algorithm parameters. This was described above.
	If `NULL`, then the default parameters are used.
=workspace=
	Points to a `Workspace` structure if preallocation is used.
	If `NULL`, then it is internally allocated and freed.
=Adata, Bdata=
	User data that is passed in to each call of `Aop` and `Bop`. For sparse matrices, typically this will be a pointer to the sparse matrix structure itself.
=SortData=
	User data that is passed in to each call of `sorter`.




---
## ShiftInvert<a name="IRA_ShiftInvert" />

This is the driver routine for shift-and-invert mode problems.
This is the driver to use for computing eigenvalues near a specified target (nearest the shift).
Both generalized and non-generalized problems are handled here, however the operator callback functions must perform different actions in the two cases, so care must be taken in making sure the action matches the problem being solved.

Solves an eigensystem `A * X = lambda * B * X`, where only matrix-vector multiplication routines are available for `A` and `B` (typically the case for large sparse matrices), and a small number of eigenvalues/vectors are needed.
The case of `B = I` is handled specially.

Note that internally, the actual problem being solved is
    inv(A - sigma * B) * B * x = nu * x
where `nu = 1/(lambda - sigma)`. Therefore for eigenvalues near `sigma`, `nu` is very large.
Hence, the sort function should almost always be `LargestMagnitude`.

Returns the number of converged eigenpairs.

### Prototype

	int ShiftInvert(
		size_t n, const std::complex<double> &shift, Ashifted Aop, Bfunc Bop,
		size_t n_wanted, size_t n_arnoldi,
		SortFunc sorter,
		std::complex<double> *w,
		std::complex<double> *v, size_t ldv,
		Params *params = NULL,
		Workspace *workspace = NULL,
		void *Adata = NULL, void *Bdata = NULL, void *SortData = NULL);

### Arguments:

=n=
	The number of rows/columns of matrices `A` and `B`.
=shift=
	The value near which the desired eigenvalues are located.
=Aop, Bop=
	If `Bop` is `NULL`, then the non-generalized problem is solved, and `Aop` must perform 
		y <- inv(A - shift*I)*x
	If `Bop` is not `NULL`, then the generalized problem is solved, and `Aop` must perform
		y <- inv(A - shift*B)*x
	and `Bop` must perform `y <- B*x`.
	The third argument is user data, and can be passed in via `Adata` and `Bdata`.
	Note that the prototype for `Aop` differs between driver routines.
=n_wanted=
	The number of desired eigenvalues/vectors to be computed (`n_wanted <= n-2`).
=n_arnoldi=
	The maximum number of Arnoldi vectors to use. It is required that `n_wanted+1 <= n_arnoldi <= n`.
	A reasonable choice if unsure is `2*n_wanted+1`.
=sorter=
	The eigenvalue sort function used for selecting which part of the spectrum to compute.
	For targeted problems, the predefined function `LargestMagnitude` should be used.
=w=
	Vector of length `n_wanted` to store the eigenvalues.
=v=
	Pointer to the first element of the matrix which will hold the eigenvectors.
	It must have at least `n` rows and `n_wanted` columns. The eigenvectors are stored in the columns in the same order as the eigenvalues.
=ldv=
	Leading dimension of matrix `v` (typically `n` unless `v` is a submatrix of a larger matrix).
=params=
	Specifies the IRA algorithm parameters. This was described above.
	If `NULL`, then the default parameters are used.
=workspace=
	Points to a `Workspace` structure if preallocation is used.
	If `NULL`, then it is internally allocated and freed.
=Adata, Bdata=
	User data that is passed in to each call of `Aop` and `Bop`. For sparse matrices, typically this will be a pointer to the sparse matrix structure itself.
=SortData=
	User data that is passed in to each call of `sorter`.

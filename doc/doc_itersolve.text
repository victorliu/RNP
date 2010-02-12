% RNP Documentation -- Iterative linear solvers
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Iterative linear solvers

This module contains several iterative linear solvers (typically used for sparse matrices).
They all share a common interface.
All functions are in namespace `RNP`.

## Contents

* [Iterative linear solver interface](#RNP_IterativeLinSolveInterface)
* [GMRES](#RNP_GMRES)
* [IDRs](#RNP_IDRs)
* [BiCGSTABl](#RNP_BiCGSTABl)

---
## Iterative linear solver interface<a name="RNP_IterativeLinSolveInterface" />

All iterative linear solvers share a common interface.
This makes it easy to swap solvers or use them interchangably as subunits of higher level algorithms.

### Prototype

	struct LinearSolverWorkspace;
	struct LinearSolverParams;
	struct NAME_params;
	int NAME(
		const size_t n,
		void (*A)(const std::complex<double>*, std::complex<double>*, void*),
		std::complex<double>* b,
		std::complex<double>* x,
		LinearSolverParams &ls_params,
		const NAME_params &params,
		LinearSolverWorkspace *workspace = NULL,
		void *Adata = NULL,
		void (*P)(std::complex<double>*, void*) = NULL,
		void *Pdata = NULL
	);
	void NAME_alloc(size_t n, const NAME_params &params, LinearSolverWorkspace *workspace);
	void NAME_free(size_t n, const NAME_params &params, LinearSolverWorkspace *workspace);

Associated with each solver are specialized structures for selecting the solver parameters, as well as specific functions for managing workspaces.

### Arguments:

=n=
	Dimension of the problem; size of matrix `A`.
=A=
	A function which multiplies the matrix `A` by the first argument
	and returns the result in the second. The second argument must
	be manually zeroed. The third parameter is user data, passed
	in through Adata.
=b=
	RHS of the problem (assumed to be length n). Possibly overwritten on exit.
=x=
	The solution vector is returned here (assumed to be length n).
=ls_params=
	A parameter structure containing options common to all linear solvers.
	On exit, each field is overwritten with the actual details of the solve.
	More information about this is described below.
=params=
	A parameter structure containing options specific to each solver.
=workspace=
	The temporary workspace for the solver. Use the provided `NAME_alloc` and `NAME_free` functions to allocate the proper size space. This parameter can be NULL, in which case the solver automatically allocates and frees the space internally.
=Adata=
	User data for function `A`.
=P=
	A preconditioner which solves `P * y = x` for `y`, overwriting `x` with `y`, where `x` is the first argument. The second argument is user data which is passed in through Pdata.
=Pdata=
	User data for function `B`.

### LinearSolverParams

A description of each field of this structure is given below.

=max_iterations=
	The maximum number of iterations allowed in the solve.
	Set to 0 to let the routine automatically select a reasonable default based on `n` and other parameters.
	On exit from the solver, this field is overwritten by the actual number of iterations required.
	The default value is 0.
=tol=
	The tolerance of the residual norm required for termination.
	The residual norm must be `tol` times the norm of the RHS for successful termination of the solve.
	On exit, this field is overwritten by the actual ratio of the norm of the residual to the norm of the RHS.
	The default value is 1e-8;
=max_mult=
	The maximum number of matrix-vector multiplies allowed.
	Specify 0 for unlimited.
	On exit, this field is overwritten by the actual number of matrix-vector multiplies performed.
	The default value is 0.
=x_initialized=
	If true, the vector `x` is assumed to be a starting guess for the solve.
	If false, the vector `x` is overwritten with zeros.
	The default value is false.



---
## GMRES<a name="RNP_GMRES" />

This is an implementation of the Generalized Minimal Residual method.
This method solves for the minimal residual in the Krylov subspace at each step of the solve.
In doing so, it requires a substantial amount of storage space (this algorithm favors fast convergence in exchange for using a huge amount of memory).

### GMRES_params

=m=
	The maximum degree of minimum residual polynomial.
	The required workspace for GMRES is proportional to `m` (approximately `n*m`).
	The default value is 100.


---
## IDRs<a name="RNP_IDRs" />

This is an implementation of the Induced Dimension Reduction method IDR(s).
When `s=1`, this method is roughly equivalent to Bi-CGSTAB.
Increasing the parameter `s` trades off larger storage space for faster convergence.

### IDRs_params

=s=
	The dimension of the shadow space.
	The required workspace for IDRs is proportional to `s` (approximately `3*n*s`).
	The default value is 4.
=angle=
	Minimum angle between subsequent Krylov vectors.
	The default value is 0.7.


---
## BiCGSTABl<a name="RNP_BiCGSTABl" />

This is an implementation of the Bi-conjugate Gradients Stabilized (l) method.
When `l=1`, this method is equivalent to Bi-CGSTAB.
Increasing the parameter `l` trades off larger storage space for faster convergence.

### BiCGSTABl_params

=l=
	Maximum degree of the Minimum Residual polynomial.
	The required workspace for BiCGSTABl is roughly proportional to `l` (goes as `3+2*(l+1)`).
	The default value is 4.
=delta=
	Fraction of initial residual norm for which the reliable update strategy takes effect.
	The default value is 1e-2.
=angle=
	Minimum angle between subsequent vectors.
	The default value is 0.7.



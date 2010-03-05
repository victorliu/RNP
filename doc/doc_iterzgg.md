% RNP Documentation -- Iterative generalized eigensystems
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# Iterative generalized eigensystems

This module contains functions related to the iterative computation of parts of the spectrum of generalized eigensystems.
All functions are in namespace `RNP`.

## Contents

* [JDQZ](#RNP_JDQZ)

---
## JDQZ<a name="RNP_JDQZ" />

Implements the Jacobi-Davidson QZ iteration method.

Solves a generalized eigensystem `A * X = lambda * B * X`, where only matrix-vector multiplication routines are available for `A` and `B` (typically the case for large sparse matrices), and a small number of eigenvalues/vectors are needed.
The eigenvalues `lambda` are given as the ratio `alpha/beta`, just like in the dense case.
This method can solve for any portion of the spectrum; all that is required is a comparator function for eigenvalues.

The JDQZ method relies on the dense QZ iteration to find small-scale Schur decompositions, as well as an iterative linear solver for solving the correction equation to compute new search directions.
Functions and adapters are provided for using the iterative solvers in RNP.

Preconditioning is optional but highly recommended as it vastly improves the convergence speed.

Returns the number of converged eigenpairs.

### Prototype

	int JDQZ(
		size_t n, size_t kmax,
		void (*Aop)(const std::complex<double> *, std::complex<double> *, void*),
		void (*Bop)(const std::complex<double> *, std::complex<double> *, void*),
		std::complex<double> *alpha, std::complex<double> *beta,
		std::complex<double> *eivec, size_t ldeivec,
		const std::complex<double> &target,
		bool (*EvalSorter)(const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, void*),
		void (*preconditioner)(std::complex<double> *, void *),
		const JDQZParams &params,
		JDQZWorkspace *workspace = NULL,
		void *Aop_data = NULL,
		void *Bop_data = NULL,
		void *EvalSorter_data = NULL,
		void *preconditioner_data = NULL
	);

### Arguments:

=n=
	The number of rows/columns of matrices `A` and `B`.
=kmax=
	The number of desired eigenvalues/vectors to be computed (`kmax <= n`).
=Aop, Bop=
	Functions which perform matrix multiplication of `A` and `B` by a vector specified in the first argument, and places it in the second argument.
	The third argument is user data, and can be passed in via `Aop_data` and `Bop_data`.
=alpha, beta=
	Vectors of length `kmax`.
	On exit, `alpha[j]/beta[j]` will be the generalized eigenvalues.

	Note: the quotients `alpha[j]/beta[j]` may easily over- or
	underflow, and `beta[j]` may even be zero.  Thus, the user
	should avoid naively computing the ratio `alpha/beta`.
	However, alpha will be always less than and usually
	comparable with `norm(A)` in magnitude, and `beta` always less
	than and usually comparable with `norm(B)`.
=eivec=
	Pointer to the first element of the matrix which will hold the eigenvectors.
	It must have at least `n` rows and `kmax` columns. The eigenvectors are stored in the columns in the same order as the eigenvalues.
=ldeivec=
	Leading dimension of matrix `eivec` (typically `n` unless `eivec` is a submatrix of a larger matrix).
=target=
	For solves which aim to compute all eigenvalues nearest to a certain number, specify the target value here.
	For non-targeted solves, it is acceptable to take `target = 0`.
	Note that this value is not passed to the preconditioner, so the preconditioner must either not require this value or obtain it through its user data paramter.
=EvalSorter=
	A comparison (less-than) function for eigenvalues. The sort order is defined such that if eigenvalue `p` is more desirable than `q`, then `EvalSorter(p,q) = true`.
	The sort function is not passed the candidate eigenvalues directly; it is passed the ratios `alpha/beta` for each eigenvalue.
	Therefore, the first two parameters correspond to `alpha` and `beta`, respectively, of the first eigenvalue candidate, and the third and fourth arguments correspond to `alpha` and `beta` of the second eigenvalue candidate.
	The last argument is user data and needed for specifying information such as the target eigenvalue for generic sort functions.
	
	Several predefined sort functions are provided for convenience, and some of them specify what `EvalSorter_data` must be.
=preconditioner=
	A function which preconditions the matrix pencil `A - target * B`.
	It is given a vector in the first argument, to which it should apply an approximation of `inv(A - target * B)`.
	For non-targeted solves, it is acceptable to take `target = 0`.
	The second argument is user data and passed in through preconditioner_data.
=params=
	Specifies the JDQZ algorithm parameters. More details about this are provided below.
=workspace=
	Pointer to a workspace structure. If ignored (by passing in `NULL`), the workspace is automatically allocated and freed internally.
	If repeated solves will be performed, it is best to pre-allocate the space and use it repeatedly.
	Use the functions `JDQZ_alloc` and `JDQZ_free` for this purpose.
=Aop_data, Bop_data=
	User data that is passed in to each call of `Aop` and `Bop`. For sparse matrices, typically this will be a pointer to the sparse matrix structure itself.
=EvalSorter_data=
	User data that is passed in to each call of `EvalSorter`. For targeted solves, typically this will be a pointer to the target value.
=preconditioner_data=
	User data that is passed in to each call of `preconditioner`. For targeted solves, typically this will be a pointer to the target value.

### Predefined sort functions

* JDQZSort_MaxReal - Solves for eigenvalues with the greatest real part.
* JDQZSort_MinReal - Solves for eigenvalues with the least real part.
* JDQZSort_MaxAbs - Solves for eigenvalues with the greatest magnitude.
* JDQZSort_MinAbs - Solves for eigenvalues with the least magnitude.
* JDQZSort_MaxImag - Solves for eigenvalues with the greatest imaginary part.
* JDQZSort_MinImag - Solves for eigenvalues with the least imaginary part.
* JDQZSort_Nearest - Solves for eigenvalues nearest some target value. The target value must be passed in as user data in the form of a pointer to a `std::complex<double>` whose value is the target.
* JDQZSort_NearestAbs - Solves for eigenvalues whose magnitude is nearest some target value. The target value must be passed in as user data in the form of a pointer to a `double` whose value is the target.

### JDQZParams

The fields in the parameter structure are described below.

=eps=
	The tolerance of the eigensolutions: `|beta * A * x - alpha * B * x| / (|alpha/beta|) < eps`.
	The default value is 1e-9.
=lock=
	Specifies the tracking parameter.
	When eigenvalues are near convergence, the tracking parameter controls how the target adaptively changes.
	The default value is 1e-9.
=min_search_space, max_search_space=
	Specifies the bounds on the dimension of the search space.
	This parameter controls the amount of temporary storage required by the algorithm.
	If 0 is specified for either bound, a default value is computed based on the desired `kmax`.
	The default values are both 0.
=max_iters=
	The maximum number of iterations allowed in the JDQZ algorithm.
=max_mult=
	The maximum number of matrix-vector multiplications allowed in the linear solver.
=testspace=
	Specifies how the testspace is generated. This parameter is very important and its selection depends on the type of solve to be performed.
	
	* = 1: Standard Petrov (`conj(alpha) * A * v + conj(beta) * B * v`)
	* = 2: Standard "variable" Petrov
	* = 3: Harmonic Petrov (`beta * A * v - alpha * B * v`)
	* = 4: The test space is the same as the search space
	* = 5: `B * v`
	* = 6: `A * v`
	
	When a targeted solve is desired, only values 1-3 should be used, with 3 being the most preferable.
	The value 4 is not implemented, and should never be used.
	When an extremal part of the spectrum is desired, 5 should be used. Other values should be avoided unless you know what you are doing.
=linear_solver=
	A function which performs a linear solve of a given matrix operator and preconditioner.
	Its prototype is as follows:

		void linear_solver(
			size_t n,
			// These next two arguments should be passed unmodified as the linear
			// solver's matrix-vector operator data. Failure to do so results in
			// undefined behavior.
			void (*op)(const std::complex<double> *, std::complex<double> *, void*),
			void *op_data,
			std::complex<double> *x, // Place solution vector here
			std::complex<double> *b, // RHS given here, can be modified
			double tol,
			size_t &max_mult,
			// These next two arguments should be passed unmodified as the linear
			// solver's preconditioner data. Failure to do so results in undefined
			// behavior.
			void (*jdqz_precon)(std::complex<double> *x, void *precon_data),
			void *jdqz_precon_data,
			void *data); // data is passed in through linear_solver_data
	
	The parameters `op`, `op_data`, `jdqz_precon`, and `jdqz_precon_data` should be passed on to the linear solver as the matrix-vector multiplication routine, and its preconditioner (and associated user data).
	As usual, `n` is the dimension of the matrix, `tol` is the solution tolerance, and `max_mult` is the maximum number of matrix-vector multiplies allowed.
	The result should be stored in `x`, while the right hand side is given in `b` and is allowed to be modified.
	The last argument `data` is user data that is passed in through `linear_solver_data` (the next argument described below).
=linear_solver_data=
	This is user data passed to `linear_solver` each time it is called. Typically, this will be a structure specifying algorithmic parameters for the linear solver.

### Linear solver adapters

Adapters are provided for the iterative solvers in RNP.
All of them require passing in a pointer to a `JDQZ_LinearSolver_adapter_data` structure in the `linear_solver_data` field of the `JDQZParams` structure.
The adapters use the default algorithmic parameters of the solvers.
The adapters are list below:

* `JDQZ_IDRs_adapter` - adapter for `IDRs`.
* `JDQZ_GMRES_adapter` - adapter for `GMRES`.

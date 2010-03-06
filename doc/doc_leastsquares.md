% RNP Documentation -- LeastSquares
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# LeastSquares

Currently this module only provides an iterative least squares solver for large sparse systems.

## Contents

* [LSQR](#RNP_LSQR)

---
## LSQR<a name="RNP_LSQR" />

LSQR  finds a solution x to the following problems:

1. Unsymmetric equations -- solve
       A*x = b

2. Linear least squares -- solve in the least-squares sense
       A*x = b

3. Damped least squares -- solve in the least-squares sense
       (   A    )*x = ( b )
       ( damp*I )     ( 0 )

where `A` is a matrix with m rows and n columns, `b` is an m-vector, and `damp` is a scalar.
The matrix `A` is intended to be large and sparse.

Note:  LSQR uses an iterative method to approximate the solution.
The number of iterations required to reach a certain accuracy
depends strongly on the scaling of the problem.  Poor scaling of
the rows or columns of `A` should therefore be avoided where
possible.

For example, in problem 1 the solution is unaltered by
row-scaling.  If a row of `A` is very small or large compared to
the other rows of `A`, the corresponding row of `( A  b )` should be
scaled up or down.

In problems 1 and 2, the solution `x` is easily recovered
following column-scaling.  Unless better information is known,
the nonzero columns of `A` should be scaled so that they all have
the same Euclidean norm (e.g., 1.0).

In problem 3, there is no freedom to re-scale if damp is
nonzero.  However, the value of `damp` should be assigned only
after attention has been paid to the scaling of `A`.

The parameter `damp` is intended to help regularize
ill-conditioned systems, by preventing the true solution from
being very large.  Another aid to regularization is provided by
the parameter acond, which may be used to terminate iterations
before the computed solution becomes very large.

Note that `x` is not an input parameter.
If some initial estimate `x0` is known and if `damp` = 0,
one could proceed as follows:

  1. Compute a residual vector     `r0 = b - A*x0`.
  2. Use LSQR to solve the system  `A*dx = r0`.
  3. Add the correction `dx` to obtain a final solution `x = x0 + dx`.

This requires that `x0` be available before and after the call
to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
to solve `A*x = b` and k2 iterations to solve `A*dx = r0`.
If `x0` is "good", `norm(r0)` will be smaller than `norm(b)`.
If the same stopping tolerances `atol` and `btol` are used for each
system, k1 and k2 will be similar, but the final solution `x0 + dx`
should be more accurate.  The only way to reduce the total work
is to use a larger stopping tolerance for the second system.
If some value `btol` is suitable for `A*x = b`, the larger value
`btol*norm(b)/norm(r0)` should be suitable for `A*dx = r0`.

Preconditioning is another way to reduce the number of iterations.
If it is possible to solve a related system `M*x = b` efficiently,
where `M` approximates `A` in some helpful way
(e.g. `M - A` has low rank or its elements are small relative to
those of `A`), LSQR may converge more rapidly on the system `A*inv(M)*z = b`,
after which `x` can be recovered by solving `M*x = z`.

NOTE: If `A` is symmetric, LSQR should not be used!
Alternatives are the symmetric conjugate-gradient method (cg)
and/or SYMMLQ.
SYMMLQ is an implementation of symmetric cg that applies to
any symmetric `A` and will converge more rapidly than LSQR.
If `A` is positive definite, there are other implementations of
symmetric cg that require slightly less work per iteration
than SYMMLQ (but will take the same number of iterations).


The following quantities are used in discussing the subroutine
parameters:

    Abar   =  (   A    ),        bbar  =  ( b )
              ( damp*I )                  ( 0 )

    r      =  b - A*x,           rbar  =  bbar - Abar*x

    rnorm  =  sqrt( norm(r)^2  +  damp^2 * norm(x)^2 )
           =  norm( rbar )

    eps    =  the relative precision of floating-point arithmetic
              on the machine being used.  On most machines,
              relpr is about 1.0e-7 and 1.0d-16 in single and double
              precision respectively.

LSQR  minimizes the function `rnorm` with respect to `x`.

### Prototype

	template <class T>
	struct LSQR{
		LSQR(
			size_t m,
			size_t n,
			void (*Aop)(bool trans, T *x, T *y, void *data),
			T damp,
			T *b,
			T *x,
			Params &params,
			T *work = NULL, // length 2*n
			void *Aop_data = NULL,
			OptionalOutputs *out = NULL
		);
	};

### Arguments:

=m=
	The number of rows of matrix `A`.
=n=
	The number of columns of matrix `A`.
=Aop=
	A function which multiples matrix `A` by a vector.
	If trans = false, compute `y = y + A*x`.
	If trans = true,  compute `x = x + A'*y`. (Conjugate transpose multiply)
	The vectors x and y are input parameters in both cases.
	If trans = true, y should be altered without changing x.
	If trans = false, x should be altered without changing y.
=damp=
	The value of the damping parameter. This is typically 0 for ordinary least squares.
	For complex numbers, the value `abs(damp)` is used.
=b=
	Pointer to the right hand side vector `b`. It is overwritten.
=x=
	Pointer to the returned (approximate) solution vector.
=params=
	A parameter structure specifying the tolerances and iteration limits of the algorithm.
	This is described in more detail below.
=work=
	Workspace used for temporary storage. If left `NULL`, then it will be internally allocated and freed.
	If non-`NULL`, must be of length `2*n`.
=Aop_data=
	Pointer to user data that is passed to `Aop`.
=out=
	Pointer to a structure holding optional output information.
	This is described in more detail below

### `Params` structure

The default constructor for the `Params` structure provides reasonable default values.
A description of each field follows.

=atol=
	An estimate of the relative error in the data defining the matrix `A`.
	For example, if `A` is accurate to about 6 digits, set `atol = 1.0e-6`.
	The default is 0, which sets `atol = eps`.
=btol=
	An estimate of the relative error in the data defining the RHS vector `b`.
	For example, if `b` is accurate to about 6 digits, set `btol = 1.0e-6`.
	The default is 0, which sets `btol = eps`.
=max_Acond=
	An upper limit on `cond(Abar)`, the apparent
	condition number of the matrix `Abar`.
	Iterations will be terminated if a computed
	estimate of `cond(Abar)` exceeds max_Acond.
	This is intended to prevent certain small or
	zero singular values of `A` or `Abar` from
	coming into effect and causing unwanted growth
	in the computed solution.
	
	`max_Acond` and `damp` may be used separately or
	together to regularize ill-conditioned systems.
	
	Normally, `max_Acond` should be in the range
	1000 to 1/`eps`.
	Suggested value:
	`max_Acond = 1/(100*eps)`  for compatible systems,
	`max_Acond = 1/(10*sqrt(eps))` for least squares.
	`eps` is the machine precision.
	
	The default is 0, which sets `max_Acond = 1/eps`.
=max_iterations=
	The maximum number of iterations allowed.
	Suggested value: `max_iterations = n/2` for well-conditioned systems with clustered singular values, `max_iterations = 4*n` otherwise.
	The default is the latter.

### `OptionalOutputs` structure

These are outputs which provide status information as well as various condition and norm estimates.

=info=
	An integer giving the reason for termination:
	=0=
		`x = 0` is the exact solution.
		No iterations were performed.
	=1=
		The equations `A*x = b` are probably compatible.
		`norm(A*x - b)` is sufficiently small, given the values of `atol` and `btol`.
	=2=
		`damp` is zero. The system `A*x = b` is probably not compatible.
		A least-squares solution has been obtained that is sufficiently accurate, given the value of `atol`.
	=3=
		`damp` is nonzero. A damped least-squares solution has been obtained that is sufficiently accurate, given the value of `atol`.
	=4=
		An estimate of `cond(Abar)` has exceeded `max_Acond`. 
		The system `A*x = b` appears to be ill-conditioned.
		Otherwise, there could be an error in the function `Aop`.
	=5=
		The iteration limit `max_iterations` was reached.
=standard_error=
	If non NULL, the dimension of `standard_error` must be `n` or more. `standard_error` then returns standard error estimates for the components of x.
	For each `i`, `standard_error[i]` is set to the value `rnorm * sqrt( sigma(i,i) / t )`, where `sigma(i,i)` is an estimate of the `i`-th diagonal of the inverse of `Abar'*Abar`.
	And:
	    t = 1      if  m <= n,
	    t = m - n  if  m > n  and  damp = 0,
	    t = m      if  damp = 0.
	
	If `NULL`, `standard_error` will not be referenced.
=anorm=
	An estimate of the Frobenius norm of `Abar`. This is the square-root of the sum of squares
	of the elements of `Abar`. If `damp` is small and if the columns of `A`
	have all been scaled to have length 1.0, `anorm` should increase to roughly `sqrt(n)`.
	A radically different value for anorm may indicate an error in the function `Aop` (there
	may be an inconsistency between the conjugate transpose and the non-transposed versions).
=acond=
	An estimate of `cond(Abar)`, the condition number of `Abar`.  A very high value of acond may again indicate an error in aprod.
=rnorm=
	An estimate of the final value of `norm(rbar)`, the function being minimized (see notation above).  This will be small if `A*x = b` has a solution.
=arnorm=
	An estimate of the final value of `norm( Abar'*rbar )`, the norm of the residual for the usual normal equations. This should be small in all cases.
	(`arnorm` will often be smaller than the true value computed from the output vector `x`.)
=xnorm=
	An estimate of the norm of the final solution vector `x`.
#ifndef _RNP_LEAST_SQUARES_H_
#define _RNP_LEAST_SQUARES_H_

#include "TBLAS.h"
#include "TLASupport.h"

namespace RNP{

template <class T>
struct LSQR{

//  LSQR  finds a solution x to the following problems:
//
//  1. Unsymmetric equations --    solve  A*x = b
//
//  2. Linear least squares  --    solve  A*x = b
//                                 in the least-squares sense
//
//  3. Damped least squares  --    solve  (   A    )*x = ( b )
//                                        ( damp*I )     ( 0 )
//                                 in the least-squares sense
//
//  where A is a matrix with m rows and n columns, b is an
//  m-vector, and damp is a scalar.
//  The matrix A is intended to be large and sparse.  It is accessed
//  by means of subroutine calls of the form
//
//             call Aop( trans, x, y)
//
//  which must perform the following functions:
//
//             If trans = false, compute  y = y + A*x.
//             If trans = true,  compute  x = x + A'*y.
//
//  The vectors x and y are input parameters in both cases.
//  If  trans = true,  y should be altered without changing x.
//  If  trans = false, x should be altered without changing y.
//  The parameters leniw, lenrw, iw, rw may be used for workspace
//  as described below. For the case of complex A, if trans = true,
//  the conjugate transpose should be applied, not the transpose.
//
//  The rhs vector b is input via u, and subsequently overwritten.
//
//
//  Note:  LSQR uses an iterative method to approximate the solution.
//  The number of iterations required to reach a certain accuracy
//  depends strongly on the scaling of the problem.  Poor scaling of
//  the rows or columns of A should therefore be avoided where
//  possible.
//
//  For example, in problem 1 the solution is unaltered by
//  row-scaling.  If a row of A is very small or large compared to
//  the other rows of A, the corresponding row of ( A  b ) should be
//  scaled up or down.
//
//  In problems 1 and 2, the solution x is easily recovered
//  following column-scaling.  Unless better information is known,
//  the nonzero columns of A should be scaled so that they all have
//  the same Euclidean norm (e.g., 1.0).
//
//  In problem 3, there is no freedom to re-scale if damp is
//  nonzero.  However, the value of damp should be assigned only
//  after attention has been paid to the scaling of A.
//
//  The parameter damp is intended to help regularize
//  ill-conditioned systems, by preventing the true solution from
//  being very large.  Another aid to regularization is provided by
//  the parameter acond, which may be used to terminate iterations
//  before the computed solution becomes very large.
//
//  Note that x is not an input parameter.
//  If some initial estimate x0 is known and if damp = 0,
//  one could proceed as follows:
//
//    1. Compute a residual vector     r0 = b - A*x0.
//    2. Use LSQR to solve the system  A*dx = r0.
//    3. Add the correction dx to obtain a final solution x = x0 + dx.
//
//  This requires that x0 be available before and after the call
//  to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
//  to solve A*x = b and k2 iterations to solve A*dx = r0.
//  If x0 is "good", norm(r0) will be smaller than norm(b).
//  If the same stopping tolerances atol and btol are used for each
//  system, k1 and k2 will be similar, but the final solution x0 + dx
//  should be more accurate.  The only way to reduce the total work
//  is to use a larger stopping tolerance for the second system.
//  If some value btol is suitable for A*x = b, the larger value
//  btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
//
//  Preconditioning is another way to reduce the number of iterations.
//  If it is possible to solve a related system M*x = b efficiently,
//  where M approximates A in some helpful way
//  (e.g. M - A has low rank or its elements are small relative to
//  those of A), LSQR may converge more rapidly on the system
//        A*M'*z = b,
//  after which x can be recovered by solving M*x = z.
//
//  NOTE: If A is symmetric, LSQR should not be used!
//  Alternatives are the symmetric conjugate-gradient method (cg)
//  and/or SYMMLQ.
//  SYMMLQ is an implementation of symmetric cg that applies to
//  any symmetric A and will converge more rapidly than LSQR.
//  If A is positive definite, there are other implementations of
//  symmetric cg that require slightly less work per iteration
//  than SYMMLQ (but will take the same number of iterations).
//
//
//  Notation
//  --------
//
//  The following quantities are used in discussing the subroutine
//  parameters:
//
//  Abar   =  (   A    ),          bbar  =  ( b )
//            ( damp*I )                    ( 0 )
//
//  r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
//
//  rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
//         =  norm( rbar )
//
//  relpr  =  the relative precision of floating-point arithmetic
//  (or eps)  on the machine being used.  On most machines,
//            relpr is about 1.0e-7 and 1.0d-16 in single and double
//            precision respectively.
//
//  LSQR  minimizes the function rnorm with respect to x.
//
//
//  Parameters
//  ----------
//
//  m       input      m, the number of rows in A.
//
//  n       input      n, the number of columns in A.
//
//  Aop     external   See above.
//
//  damp    input      The damping parameter for problem 3 above.
//                     (damp should be 0.0 for problems 1 and 2.)
//                     If the system A*x = b is incompatible, values
//                     of damp in the range 0 to sqrt(relpr)*norm(A)
//                     will probably have a negligible effect.
//                     Larger values of damp will tend to decrease
//                     the norm of x and reduce the number of 
//                     iterations required by LSQR.
//
//                     The work per iteration and the storage needed
//                     by LSQR are the same for all values of damp.
//
//  u(m)    input      The rhs vector b.  Beware that u is
//                     over-written by LSQR.
//
//  x(n)    output     Returns the computed solution x.
//
//
////
//  References
//  ----------
//
//  C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
//       linear equations and sparse least squares,
//       ACM Transactions on Mathematical Software 8, 1 (March 1982),
//       pp. 43-71.
//
//  C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
//       linear equations and least-squares problems,
//       ACM Transactions on Mathematical Software 8, 2 (June 1982),
//       pp. 195-209.
//
//  C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
//       Basic linear algebra subprograms for Fortran usage,
//       ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
//       pp. 308-323 and 324-325.
//  ------------------------------------------------------------------
//

struct Params{
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	
	// An estimate of the relative error in the data
	// defining the matrix A.  For example,
	// if A is accurate to about 6 digits, set
	// atol = 1.0e-6.
	real_type atol;

	// An estimate of the relative error in the data
	// defining the rhs vector b.  For example,
	// if b is accurate to about 6 digits, set
	// btol = 1.0e-6 .
	real_type btol;
	
	// An upper limit on cond(Abar), the apparent
	// condition number of the matrix Abar.
	// Iterations will be terminated if a computed
	// estimate of cond(Abar) exceeds max_Acond.
	// This is intended to prevent certain small or
	// zero singular values of A or Abar from
	// coming into effect and causing unwanted growth
	// in the computed solution.
	//
	// max_Acond and damp may be used separately or
	// together to regularize ill-conditioned systems.
	//
	// Normally, max_Acond should be in the range
	// 1000 to 1/eps.
	// Suggested value:
	// max_Acond = 1/(100*eps)  for compatible systems,
	// max_Acond = 1/(10*sqrt(eps)) for least squares.
	real_type max_Acond;
	
	// Note:  If the user is not concerned about the parameters
	// atol, btol and max_Acond, any or all of them may be set
	// to zero.  The effect will be the same as the values
	// eps, eps, and 1/eps respectively. (eps is machine precision)

	size_t max_iterations; // max iterations, on exit, actual iterations used
	Params(){
		// 0 for these three parameters indicates machine precision should be used
		atol = 0;
		btol = 0;
		max_Acond = 0;
		
		max_iterations = 0; // autoselect, use n/2 if the system is "nice", otherwise use 4*n
	}
};

struct OptionalOutputs{
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	
	// info: an integer giving the reason for termination:
	//
	// 0       x = 0  is the exact solution.
	//         No iterations were performed.
	//
	// 1       The equations A*x = b are probably
	//         compatible.  Norm(A*x - b) is sufficiently
	//         small, given the values of atol and btol.
	//
	// 2       damp is zero.  The system A*x = b is probably
	//         not compatible.  A least-squares solution has
	//         been obtained that is sufficiently accurate,
	//         given the value of atol.
	//
	// 3       damp is nonzero.  A damped least-squares
	//         solution has been obtained that is sufficiently
	//         accurate, given the value of atol.
	//
	// 4       An estimate of cond(Abar) has exceeded
	//         max_Acond.  The system A*x = b appears to be
	//         ill-conditioned.  Otherwise, there could be an
	//         error in subroutine aprod.
	//
	// 5       The iteration limit max_iterations was reached.
	int info;
	
	// If non NULL, the dimension of standard_error must be n or more. standard_error then returns standard error
	// estimates for the components of x. For each i, standard_error[i] is set to the value
	//    rnorm * sqrt( sigma(i,i) / t ),
	// where sigma(i,i) is an estimate of the i-th diagonal of the inverse of Abar'*Abar
	// and  t = 1      if  m <= n,
	//      t = m - n  if  m > n  and  damp = 0,
	//      t = m      if  damp = 0.
	//
	// If NULL, standard_error will not be referenced.
	real_type *standard_error; // leave NULL if not wanted
	
	// An estimate of the Frobenius norm of  Abar. This is the square-root of the sum of squares
	// of the elements of Abar. If damp is small and if the columns of A
	// have all been scaled to have length 1.0, anorm should increase to roughly sqrt(n).
	// A radically different value for anorm may indicate an error in subroutine aprod (there
	// may be an inconsistency between modes 1 and 2).
	real_type anorm;
	
	// An estimate of cond(Abar), the condition number of Abar.  A very high value of acond
	// may again indicate an error in aprod.
	real_type acond;

	// An estimate of the final value of norm(rbar), the function being minimized (see notation
	// above).  This will be small if A*x = b has a solution.
	real_type rnorm;

	// An estimate of the final value of norm( Abar'*rbar ), the norm of
	// the residual for the usual normal equations. This should be small in all cases.
	// (arnorm will often be smaller than the true value computed from the output vector x.)
	real_type arnorm;
	
	// An estimate of the norm of the final solution vector x.
	real_type xnorm;
	
	OptionalOutputs():standard_error(NULL){}
};

LSQR(
	size_t m, // rows of A
	size_t n, // cols of A
	void (*Aop)(bool trans, T *x, T *y, void *),
	T damp,   // typically 0 for usual cases, only abs(damp) matters
	T *u,     // RHS, len = m
	T *x,     // solution output, len = n
	Params &params,
	T *work = NULL, // length 2*n
	void *Aop_data = NULL,
	OptionalOutputs *out = NULL
){
	typedef typename RNP::TBLAS::_RealOrComplexChooser<T>::real_type real_type;
	
	// Set up workspace
	T *v, *w;
	if(work == NULL){
		v = new T[2*n];
		w = v+n;
	}

	//	Local copies of output variables.  Output vars are assigned at exit.
	int	istop  = 0;
	size_t itn = 0;
	real_type
		anorm  = 0.,
		acond  = 0.,
		rnorm  = 0.,
		arnorm = 0.,
		xnorm  = 0.;
	real_type rdamp = RNP::TBLAS::_RealOrComplexChooser<T>::_abs(damp);

	//	Local variables

	const bool damped = (rdamp > 0);
	const bool wantse = (NULL != out && NULL != out->standard_error);
	real_type *se = NULL;
	if(wantse){
		se = out->standard_error;
	}
	
	size_t i, nconv;
	real_type
		alfopt, alpha, arnorm0, beta, bnorm,
		cs, cs1, cs2, ctol,
		delta, dknorm, dnorm, dxk, dxmax,
		gamma, gambar, phi, phibar, psi,
		res2, rho, rhobar, rhbar1,
		rhs, rtol, sn, sn1, sn2,
		tau, temp, test1, test2, test3,
		theta, t1, t2, t3, xnorm1, z, zbar;

	size_t nstop  =   0;
	ctol   =   0.;
	if(params.max_Acond > 0.) ctol = 1. / params.max_Acond;
	anorm  =   0.;
	acond  =   0.;
	dnorm  =   0.;
	dxmax  =   0.;
	res2   =   0.;
	psi	   =   0.;
	xnorm  =   0.;
	xnorm1 =   0.;
	cs2	   = - 1.;
	sn2	   =   0.;
	z	   =   0.;

	//	------------------------------------------------------------------
	//	Set up the first vectors u and v for the bidiagonalization.
	//	These satisfy  beta*u = b,	alpha*v = A'*u.
	//	------------------------------------------------------------------
	for(size_t i = 0; i < n; i++){ v[i] = 0; }
	for(size_t i = 0; i < n; i++){ x[i] = 0; }

	if(wantse){
		for(size_t i = 0; i < n; i++){ se[i] = 0; }
	}
	
	alpha  =   0.;
	beta   =   RNP::TBLAS::Norm2(m, u, 1);

	if(beta > 0.){
		RNP::TBLAS::Scale( m, (1. / beta), u, 1 );
		Aop(true, v, u, Aop_data);
		alpha  =   RNP::TBLAS::Norm2( n, v, 1 );
	}

	if(alpha > 0.){
		RNP::TBLAS::Scale( n, (1. / alpha), v, 1 );
		RNP::TBLAS::Copy( n, v, 1, w, 1 );
	}

	arnorm = arnorm0 = alpha * beta;
	if(arnorm != 0.){
		rhobar = alpha;
		phibar = beta;
		bnorm  = beta;
		rnorm  = beta;

		const size_t max_iter = (0 == params.max_iterations ? 4*n : params.max_iterations);
		//	==================================================================
		//	Main iteration loop.
		//	==================================================================
		while(1){
			++itn;
			
			// ------------------------------------------------------------------
			// Perform the next step of the bidiagonalization to obtain the
			// next	 beta, u, alpha, v.	 These satisfy the relations
			//			  beta*u  =	 A*v  -	 alpha*u,
			//			 alpha*v  =	 A'*u	 -	beta*v.
			// ------------------------------------------------------------------
			RNP::TBLAS::Scale( m, (- alpha), u, 1 );
			Aop(false, v, u, Aop_data);
			beta   =   RNP::TBLAS::Norm2( m, u, 1 );

			// Accumulate  anorm = || Bk ||
			//					 =	sqrt( sum of  alpha**2 + beta**2 + damp**2 ).

			temp   =   RNP::TLASupport::Pythag3( alpha, beta, rdamp );
			anorm  =   RNP::TLASupport::Pythag2( anorm, temp );

			if(beta > 0.){
				RNP::TBLAS::Scale( m, (1. / beta), u, 1 );
				RNP::TBLAS::Scale( n, (- beta), v, 1 );
				Aop(true, v, u, Aop_data);
				alpha  =   RNP::TBLAS::Norm2( n, v, 1 );
				if(alpha > 0.){
					RNP::TBLAS::Scale( n, (1. / alpha), v, 1 );
				}
			}

			// ------------------------------------------------------------------
			// Use a plane rotation to eliminate the damping parameter.
			// This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
			// ------------------------------------------------------------------
			rhbar1 = rhobar;
			if(damped){
				rhbar1 = RNP::TLASupport::Pythag2( rhobar, rdamp );
				cs1	   = rhobar / rhbar1;
				sn1	   = rdamp	/ rhbar1;
				psi	   = sn1 * phibar;
				phibar = cs1 * phibar;
			}

			// ------------------------------------------------------------------
			// Use a plane rotation to eliminate the subdiagonal element (beta)
			// of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
			// ------------------------------------------------------------------
			rho	   =   RNP::TLASupport::Pythag2( rhbar1, beta );
			cs	   =   rhbar1 / rho;
			sn	   =   beta	  / rho;
			theta  =   sn * alpha;
			rhobar = - cs * alpha;
			phi	   =   cs * phibar;
			phibar =   sn * phibar;
			tau	   =   sn * phi;

			// ------------------------------------------------------------------
			// Update  x, w	 and (perhaps) the standard error estimates.
			// ------------------------------------------------------------------
			t1 =   phi	 / rho;
			t2 = - theta / rho;
			t3 =   1.	/ rho;
			dknorm =   0.;

			if(wantse){
				for(i = 0; i < n; i++){
					T t	   =  w[i];
					x[i]   =  t1*t	+  x[i];
					w[i]   =  t2*t	+  v[i];
					real_type rt = RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(t3*t);
					se[i] += rt;
					dknorm += rt;
				}
			}else{
				for(i = 0; i < n; i++){
					T t	   =  w[i];
					x[i]   =  t1*t	+  x[i];
					w[i]   =  t2*t	+  v[i];
					dknorm += RNP::TBLAS::_RealOrComplexChooser<T>::_abs2(t3*t);;
				}
			}

			// ------------------------------------------------------------------
			// Monitor the norm of d_k, the update to x.
			// dknorm = norm( d_k )
			// dnorm  = norm( D_k ),		where	D_k = (d_1, d_2, ..., d_k )
			// dxk	  = norm( phi_k d_k ),	where new x = x_k + phi_k d_k.
			// ------------------------------------------------------------------
			dknorm = sqrt( dknorm );
			dnorm  = RNP::TLASupport::Pythag2( dnorm, dknorm );
			dxk	   = fabs( phi * dknorm );
			if(dxmax < dxk){
				dxmax	=  dxk;
			}

			// ------------------------------------------------------------------
			// Use a plane rotation on the right to eliminate the
			// super-diagonal element (theta) of the upper-bidiagonal matrix.
			// Then use the result to estimate	norm(x).
			// ------------------------------------------------------------------
			delta  =   sn2 * rho;
			gambar = - cs2 * rho;
			rhs    =   phi - delta * z;
			zbar   =   rhs / gambar;
			xnorm  =   RNP::TLASupport::Pythag2( xnorm1, zbar  );
			gamma  =   RNP::TLASupport::Pythag2( gambar, theta );
			cs2	   =   gambar / gamma;
			sn2	   =   theta / gamma;
			z      =   rhs / gamma;
			xnorm1 =   RNP::TLASupport::Pythag2( xnorm1, z	   );

			// ------------------------------------------------------------------
			// Test for convergence.
			// First, estimate the norm and condition of the matrix	 Abar,
			// and the norms of	 rbar  and	Abar'*rbar.
			// ------------------------------------------------------------------
			acond = anorm * dnorm;
			res2 = RNP::TLASupport::Pythag2( res2 , psi	   );
			rnorm = RNP::TLASupport::Pythag2( res2 , phibar );
			arnorm = alpha * fabs( tau );

			// Now use these norms to estimate certain other quantities,
			// some of which will be small near a solution.

			alfopt = sqrt( rnorm / (dnorm * xnorm) );
			test1 = rnorm /	bnorm;
			test2 = 0.;
			if(rnorm > 0.) test2 = arnorm / (anorm * rnorm);
			// if(arnorm0 > 0.) test2 = arnorm / arnorm0;  //(Michael Friedlander's modification)
			test3 = 1.	/  acond;
			t1 = test1 / (1. +  anorm * xnorm / bnorm);
			rtol = params.btol + params.atol * anorm * xnorm / bnorm;

			// The following tests guard against extremely small values of
			// atol, btol  or  ctol.  (The user may have set any or all of
			// the parameters  atol, btol, max_Acond  to zero.)
			// The effect is equivalent to the normal tests using
			// atol = relpr,  btol = relpr,	 max_Acond = 1/relpr.

			t3 = 1. + test3;
			t2 = 1. + test2;
			t1 = 1. + t1;
			if(itn >= max_iter) istop = 5;
			if(t3  <= 1.) istop = 4;
			if(t2  <= 1.) istop = 2;
			if(t1  <= 1.) istop = 1;

			// Allow for tolerances set by the user.

			if(test3 <= ctol) istop = 4;
			if(test2 <= params.atol) istop = 2;
			if(test1 <= rtol) istop = 1;   //(Michael Friedlander had this commented out)

			// ------------------------------------------------------------------
			// Stop if appropriate.
			// The convergence criteria are required to be met on  nconv
			// consecutive iterations, where  nconv	 is set below.
			// Suggested value:	 nconv = 1, 2  or  3.
			// ------------------------------------------------------------------

			if(istop == 0){
				nstop  = 0;
			}else{
				nconv  = 1;
				++nstop;
				if(nstop < nconv  &&  itn < params.max_iterations) istop = 0;
			}

			if(istop != 0) break;
			
		}
		//	==================================================================
		//	End of iteration loop.
		//	==================================================================

		//	Finish off the standard error estimates.

		if(wantse){
			real_type t = 1.;
			if(m > n){ t = m - n; }
			if(damped ){ t = m; }
			t = rnorm / sqrt(t);
		  
			for(i = 0; i < n; i++){
				se[i]  = t * sqrt(se[i]);
			}
		}

	}
	//	Decide if istop = 2 or 3.
	if(damped && istop == 2) istop = 3;

	//	Assign output variables from local copies.
	params.max_iterations = itn;
	if(NULL != out){
		out->info   = istop;
		out->anorm  = anorm;
		out->acond  = acond;
		out->rnorm  = rnorm;
		out->arnorm = test2;
		out->xnorm  = xnorm;
	}
	if(work == NULL){
		delete [] v;
	}
}

}; // namespace LSQR
}; // namespace RNP

#endif // _RNP_LEAST_SQUARES_H_

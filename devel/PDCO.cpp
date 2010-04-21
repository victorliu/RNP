// From Clp-1.11.1 of Coin-OR

#include <cstddef>
#include "Optimization.h"
#include "LeastSquares.h"
#include "TBLAS.h"
#include <cstdio>

typedef struct{
  double  atolmin;
  double  r3norm;
  double  LSdamp;
  double* deltay;
} Info;
typedef struct{
  double  atolold;
  double  atolnew;
  double  r3ratio;
  int   istop;
  int   itncg;
} Outfo;


static void Ap(bool trans, size_t m, size_t n, RNP::Optimization::PDCO::MatrixFunction A, double *diag1_, double diag2_, double *x, double *y, void *data){
  double *temp = new double[n]; RNP::TBLAS::Fill(n, 0., temp, 1);
  if (false == trans){
    A( true, m, n, temp, y, data);
    for (int k=0; k<n; k++)
      x[k] += (diag1_[k] * temp[k]);
    for (int k=0; k<m; k++)
      x[n+k] += (diag2_ * y[k]);
  }else{
    for (int k=0; k<n; k++)
      temp[k] = diag1_[k] * y[k];
    A(false, m, n, x, temp, data);
    for (int k=0; k<m; k++)
      x[k] += diag2_ * y[n+k];
  }
  delete [] temp;
  return;
}

static void do_lsqr(size_t m, size_t n, RNP::Optimization::PDCO::MatrixFunction A, double *diag1, double diag2, double *b,  
		double damp, double atol, double btol, double conlim, int itnlim, 
		bool show, Info info, double *x , int *istop,
		int *itn, Outfo *outfo, bool precon, double *Pr, void *data){

 /**
 Special version of LSQR for use with pdco.m.
 It continues with a reduced atol if a pdco-specific test isn't
 satisfied with the input atol.
*/

//     Initialize.
  static char term_msg[8][80] = 
  {
    "The exact solution is x = 0",
    "The residual Ax - b is small enough, given ATOL and BTOL",
    "The least squares error is small enough, given ATOL",
    "The estimated condition number has exceeded CONLIM",
    "The residual Ax - b is small enough, given machine precision",
    "The least squares error is small enough, given machine precision",
    "The estimated condition number has exceeded machine precision",
    "The iteration limit has been reached"
  };

  //  printf("***************** Entering LSQR *************\n");
  
  char str1[100], str2[100], str3[100], str4[100], head1[100], head2[100];

  *itn     = 0;            *istop   = 0;       int nstop  = 0;
  double ctol   = 0;         
  if (conlim > 0) ctol = 1/conlim;

  double anorm  = 0;         double acond  = 0;
  double dampsq = damp*damp; double ddnorm = 0;      double res2   = 0;
  double xnorm  = 0;         double xxnorm = 0;      double z      = 0;
  double cs2    = -1;        double sn2    = 0;

  // Set up the first vectors u and v for the bidiagonalization.
  // These satisfy  beta*u = b,  alfa*v = A'u.

  double *u = new double[m]; RNP::TBLAS::Copy(m, b, 1, u, 1);
  double *v = new double[n]; RNP::TBLAS::Fill(n, 0., v, 1);
  RNP::TBLAS::Fill(n, 0., x, 1);
  double alfa   = 0;      double beta = RNP::TBLAS::Norm2(m, u, 1);
  if (beta > 0){
    RNP::TBLAS::Scale(m, 1./beta, u, 1);
    Ap(true, m, n, A, diag1, diag2, v, u, data );
    //if (precon) v = v*Pr;
    alfa = RNP::TBLAS::Norm2(n, v, 1);
  }
  if (alfa > 0){
    RNP::TBLAS::Scale(n, 1./alfa, v, 1);
  }
  double *w = new double[n]; RNP::TBLAS::Copy(n, v, 1, u, 1);

  double arnorm = alfa * beta; 
  if (arnorm == 0){
     printf("  %s\n\n", term_msg[0]);
     return;
  }

  double rhobar = alfa;  double phibar = beta;   double bnorm  = beta;
  double rnorm  = beta;
  sprintf(head1,"   Itn      x(1)      Function");
  sprintf(head2," Compatible   LS      Norm A   Cond A");

  if (show){
    printf(" %s%s\n",head1, head2);
    double test1  = 1;  double test2  = alfa / beta;
    sprintf(str1,"%6d %12.5e %10.3e",   *itn, x[0], rnorm );
    sprintf(str2,"  %8.1e  %8.1e",       test1, test2 );
    printf("%s%s\n",str1, str2);
  }

  //----------------------------------------------------------------
  // Main iteration loop.
  //----------------------------------------------------------------
  while (*itn < itnlim){
    *itn += 1;
    // Perform the next step of the bidiagonalization to obtain the
    // next beta, u, alfa, v.  These satisfy the relations
    // beta*u  =  a*v   -  alfa*u,
    // alfa*v  =  A'*u  -  beta*v.

    RNP::TBLAS::Scale(m, -alfa, u, 1);
    //if (precon){
    //  CoinDenseVector<double> pv(v*Pr);
    //  matVecMult( 1, u, pv);
    //}else{
      Ap( false, m, m, A, diag1, diag2, u, v, data);
    //}
    beta = RNP::TBLAS::Norm2(m, u, 1);
    if (beta > 0){
      RNP::TBLAS::Scale(m, 1./beta, u, 1);
      anorm = sqrt(anorm*anorm + alfa*alfa + beta*beta +  damp*damp);
      RNP::TBLAS::Scale(n, -beta, v, 1);
      double* vv = new double[n];
      RNP::TBLAS::Fill(n, 0., vv, 1);
      Ap( true, m, n, A, diag1, diag2, vv, u, data );
      //if (precon) vv = vv*Pr;
      RNP::TBLAS::Axpy(n, 1., vv, 1, v, 1);
      delete [] vv;
      alfa  = RNP::TBLAS::Norm2(n, v, 1);
      if (alfa > 0) RNP::TBLAS::Scale(n, 1./alfa, v, 1);
    }

    // Use a plane rotation to eliminate the damping parameter.
    // This alters the diagonal (rhobar) of the lower-bidiagonal matrix.

    double rhobar1 = sqrt(rhobar*rhobar + damp*damp);
    double cs1     = rhobar / rhobar1;
    double sn1     = damp   / rhobar1;
    double psi     = sn1 * phibar;
    phibar       *= cs1;

    // Use a plane rotation to eliminate the subdiagonal element (beta)
    // of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

    double rho     = sqrt(rhobar1*rhobar1 +  beta*beta);
    double cs      =   rhobar1/ rho;
    double sn      =   beta   / rho;
    double theta   =   sn * alfa;
    rhobar  = - cs * alfa;
    double phi     =   cs * phibar;
    phibar  =   sn * phibar;
    double tau     =   sn * phi;

    // Update x and w.

    double t1      =   phi  /rho;
    double t2      = - theta/rho;
    //    dk           =   ((1/rho)*w);

    double w_norm = RNP::TBLAS::Norm2(n, w, 1);
    RNP::TBLAS::Axpy(n, t1, w, 1, x, 1);
    for(int k = 0; k < n; ++k)
		w[k]       = v[k]      +  t2*w[k];
    ddnorm  = ddnorm +  (w_norm/rho)*(w_norm/rho);
    // if wantvar, var = var  +  dk.*dk; end

    // Use a plane rotation on the right to eliminate the
    // super-diagonal element (theta) of the upper-bidiagonal matrix.
    // Then use the result to estimate  norm(x).

    double delta   =   sn2 * rho;
    double gambar  = - cs2 * rho;
    double rhs     =   phi  -  delta * z;
    double zbar    =   rhs / gambar;
    xnorm   =   sqrt(xxnorm + zbar*zbar);
    double gamma   =   sqrt(gambar*gambar + theta*theta);
    cs2     =   gambar / gamma;
    sn2     =   theta  / gamma;
    z       =   rhs    / gamma;
    xxnorm  =   xxnorm  +  z*z;

    // Test for convergence.
    // First, estimate the condition of the matrix  Abar,
    // and the norms of  rbar  and  Abar'rbar.

    acond   =   anorm * sqrt( ddnorm );
    double res1    =   phibar*phibar;
    double res2    =   res1  +  psi*psi;
    rnorm   =   sqrt( res1 + res2 );
    arnorm  =   alfa * fabs( tau );

    // Now use these norms to estimate certain other quantities,
    // some of which will be small near a solution.

    double test1   =   rnorm / bnorm;
    double test2   =   arnorm/( anorm * rnorm );
    double test3   =       1 / acond;
    t1      =   test1 / (1    +  anorm * xnorm / bnorm);
    double rtol    =   btol  +  atol *  anorm * xnorm / bnorm;

    // The following tests guard against extremely small values of
    // atol, btol  or  ctol.  (The user may have set any or all of
    // the parameters  atol, btol, conlim  to 0.)
    // The effect is equivalent to the normal tests using
    // atol = eps,  btol = eps,  conlim = 1/eps.

    if (*itn >= itnlim)
      *istop = 7; 
    if (1 + test3  <= 1)
      *istop = 6;
    if (1 + test2  <= 1)
      *istop = 5;
    if (1 + t1     <= 1)
      *istop = 4;

    // Allow for tolerances set by the user.

    if  (test3 <= ctol)
      *istop = 3;
    if  (test2 <= atol)
      *istop = 2;
    if  (test1 <= rtol)
      *istop = 1;

    //-------------------------------------------------------------------
    // SPECIAL TEST THAT DEPENDS ON pdco.m.
    // Aname in pdco   is  iw in lsqr.
    // dy              is  x
    // Other stuff     is in info.

    // We allow for diagonal preconditioning in pdDDD3.
    //-------------------------------------------------------------------
    if (*istop > 0){
      double r3new     = arnorm;
      double r3ratio   = r3new / info.r3norm;
      double atolold   = atol;
      double atolnew   = atol;
         
      if (atol > info.atolmin){
        if     (r3ratio <= 0.1){   // dy seems good
          // Relax
        }else if (r3ratio <= 0.5){ // Accept dy but make next one more accurate.
          atolnew = atolnew * 0.1;
        }else{                     // Recompute dy more accurately
	  if (show){
	    printf("\n                                ");
	    printf("                                \n");
	    printf(" %5.1f%7d%7.3f", log10(atolold), *itn, r3ratio);
	  }
          atol    = atol * 0.1;
          atolnew = atol;
          *istop   = 0;
        }
    
      outfo->atolold = atolold;
      outfo->atolnew = atolnew;
      outfo->r3ratio = r3ratio;
    }
      
    //-------------------------------------------------------------------
    // See if it is time to print something.
    //-------------------------------------------------------------------
    int prnt = 0;
    if (n     <= 40       ) prnt = 1;
    if (*itn   <= 10       ) prnt = 1;
    if (*itn   >= itnlim-10) prnt = 1;
    if (fmod(*itn,10) == 0  ) prnt = 1;
    if (test3 <=  2*ctol  ) prnt = 1;
    if (test2 <= 10*atol  ) prnt = 1;
    if (test1 <= 10*rtol  ) prnt = 1;
    if (*istop !=  0       ) prnt = 1;

    if (prnt == 1){
      if (show){
        sprintf(str1, "   %6d %12.5e %10.3e",   *itn, x[0], rnorm );
        sprintf(str2, "  %8.1e %8.1e",       test1, test2 );
        sprintf(str3, " %8.1e %8.1e",        anorm, acond );
        printf("%s%s%s\n", str1, str2, str3);
      }
    }
    if (*istop > 0)
      break;
    }
  }
  // End of iteration loop.
  // Print the stopping condition.

  if (show){
    printf("\n LSQR finished\n");
  //    disp(msg(istop+1,:))
  //    disp(' ')
    printf("%s\n", term_msg[*istop]);
    sprintf(str1, "istop  =%8d     itn    =%8d",      *istop, *itn    );
    sprintf(str2, "anorm  =%8.1e   acond  =%8.1e",  anorm, acond  );
    sprintf(str3, "rnorm  =%8.1e   arnorm =%8.1e",  rnorm, arnorm );
    sprintf(str4, "bnorm  =%8.1e   xnorm  =%8.1e",  bnorm, xnorm  );
    printf("%s %s\n", str1, str2);
    printf("%s %s\n", str3, str4);
  }
  delete [] w;
  delete [] u;
  delete [] v;
}

static inline double pdxxxmerit(size_t m, size_t n, int nlow, int nupp, int *low, int *upp, double *r1,
		double *r2, double *rL,
		double *rU, double *cL,
		double *cU ){

	// Evaluate the merit function for Newton's method.
	// It is the 2-norm of the three sets of residuals.
	double sum1, sum2;
	double f[6];
	f[0] = RNP::TBLAS::Norm2(m, r1, 1);
	f[1] = RNP::TBLAS::Norm2(n, r2, 1);
	sum1 = sum2 = 0.0;
	for (int k=0; k<nlow; k++){
		sum1 += rL[low[k]]*rL[low[k]];
		sum2 += cL[low[k]]*cL[low[k]];
	}
	f[2] = sqrt(sum1);
	f[4] = sqrt(sum2);
	sum1 = sum2 = 0.0;
	for (int k=0; k<nupp; k++){
		sum1 += rL[upp[k]]*rL[upp[k]];
		sum2 += cL[upp[k]]*cL[upp[k]];
	}
	f[3] = sqrt(sum1);
	f[5] = sqrt(sum2);

	return RNP::TBLAS::Norm2(6, f, 1);
}
inline void pdxxxresid1(size_t m, size_t n, RNP::Optimization::PDCO::MatrixFunction A, const int nlow, const int nupp, const int nfix,
		 int *low, int *upp, int *fix,
		 double *b, double *bl, double *bu, double d1, double d2,
		 double *grad, double *rL,
		 double *rU, double *x,
		 double *x1, double *x2,
		 double *y,  double *z1,
		 double *z2, double *r1,
		 double *r2, double *Pinf, double *Dinf, void *data){

	// Form residuals for the primal and dual equations.
	// rL, rU are output, but we input them as full vectors
	// initialized (permanently) with any relevant zeros.

	for (int k=0; k<nfix; k++) x[fix[k]]  = 0;

	RNP::TBLAS::Fill(m, 0., r1, 1);
	RNP::TBLAS::Fill(n, 0., r2, 1);
	A(false, m, n, r1, x, data);
	A(true, m, n, r2, y, data);
	for (int k=0; k<nfix; k++)
		r2[fix[k]]  = 0;


	for(int k = 0; k < m; ++k)
		r1[k]      = b[k]    - r1[k] - d2*d2*y[k];
	for(int k = 0; k < n; ++k)
		r2[k]      = grad[k] - r2[k] - z1[k];              // grad includes d1*d1*x
	if(nupp > 0)       
		RNP::TBLAS::Axpy(n, 1., z2, 1, r2, 1);

	for (int k=0; k<nlow; k++)
		rL[low[k]] = bl[low[k]] - x[low[k]] + x1[low[k]];
	for (int k=0; k<nupp; k++)
		rU[upp[k]] = - bu[upp[k]] + x[upp[k]] + x2[upp[k]];

	double normL = 0.0;
	double normU = 0.0;
	for (int k=0; k<nlow; k++)
		if (rL[low[k]] > normL) normL = rL[low[k]];
	for (int k=0; k<nupp; k++)
		if (rU[upp[k]] > normU) normU = rU[upp[k]];

	*Pinf    = std::max(normL, normU);  
	*Pinf    = std::max( r1[RNP::TBLAS::MaximumIndex(m, r1, 1)] , *Pinf );
	*Dinf    = r2[RNP::TBLAS::MaximumIndex(n, r2, 1)];
	*Pinf    = std::max( *Pinf, 1e-99 );
	*Dinf    = std::max( *Dinf, 1e-99 );
}

static inline void pdxxxresid2(double mu, int nlow, int nupp, int *low, int *upp,
		 double *cL, double *cU,
		 double *x1, double *x2,
		 double *z1, double *z2,
		 double *center, double *Cinf, double *Cinf0){

	// Form residuals for the complementarity equations.
	// cL, cU are output, but we input them as full vectors
	// initialized (permanently) with any relevant zeros.
	// Cinf  is the complementarity residual for X1 z1 = mu e, etc.
	// Cinf0 is the same for mu=0 (i.e., for the original problem).

	double maxXz = -1e20;
	double minXz = 1e20;

	for (int k=0; k<nlow; k++){
		double x1z1    = x1[low[k]] * z1[low[k]];
		cL[low[k]] = mu - x1z1;
		if(x1z1 > maxXz) maxXz = x1z1;
		if(x1z1 < minXz) minXz = x1z1;
	}

	for (int k=0; k<nupp; k++){
		double x2z2    = x2[upp[k]] * z2[upp[k]];
		cU[upp[k]] = mu - x2z2;
		if(x2z2 > maxXz) maxXz = x2z2;
		if(x2z2 < minXz) minXz = x2z2;
	}

	maxXz   = std::max( maxXz, 1e-99 );
	minXz   = std::max( minXz, 1e-99 );
	*center  = maxXz / minXz;

	double normL = 0.0;
	double normU = 0.0;
	for (int k=0; k<nlow; k++)
		if (cL[low[k]] > normL) normL = cL[low[k]];
	for (int k=0; k<nupp; k++)
		if (cU[upp[k]] > normU) normU = cU[upp[k]];
	*Cinf    = std::max( normL, normU);
	*Cinf0   = maxXz;
}
static inline double  pdxxxstep(size_t n, double *x, double *dx ){

	// Assumes x > 0.
	// Finds the maximum step such that x + step*dx >= 0.

	double step     = 1e+20;

	for (int k=0; k<n; k++)
		if (dx[k] < 0)
			if ((x[k]/(-dx[k])) < step)
				step = x[k]/(-dx[k]);
	return step;
}
static inline double  pdxxxstep(size_t n, int nset, int *set, double *x, double *dx ){

	// Assumes x > 0.
	// Finds the maximum step such that x + step*dx >= 0.

	double step     = 1e+20;

	for (int k=0; k<n; k++)
		if (dx[k] < 0)
			if ((x[k]/(-dx[k])) < step)
				step = x[k]/(-dx[k]);
	return step;
}
int RNP::Optimization::PDCO::Minimize(
	size_t m,
	size_t n,
	RNP::Optimization::PDCO::ConvexFunction phi,
	RNP::Optimization::PDCO::MatrixFunction A,
	double *b,  // b vector of equality constraint, length m, modified on exit
	double *bl, // lower bounds, modified on exit
	double *bu, // upper bounds, modified on exit
	const double *D1, // may be just a scalar, depending on options
	const double *D2, // may be just a scalar, depending on options
	double *x, // primal solution, length n
	double *y, // dual solution of equality constraint, length m
	double *z, // dual solution of bounds, length n
	RNP::Optimization::PDCO::Params *params,
	double *workspace,
	void *data
){
	double inf = 1.0e30;
	double eps = 1.0e-15;
	double atolold, r3ratio, Pinf, Dinf, Cinf, Cinf0;

	if(!params->x_initialized){
		RNP::TBLAS::Fill(n, 0., x, 1);
	}
	if(!params->y_initialized){
		RNP::TBLAS::Fill(m, 0., y, 1);
	}
	if(!params->z_initialized){
		RNP::TBLAS::Fill(n, 0., z, 1);
	}

	bool ifexplicit = true;

	double normb  = b[RNP::TBLAS::MaximumIndex(m, b, 1)];
	double normx0 = x[RNP::TBLAS::MaximumIndex(n, x, 1)];
	double normy0 = y[RNP::TBLAS::MaximumIndex(m, y, 1)];
	double normz0 = z[RNP::TBLAS::MaximumIndex(n, z, 1)];

	printf("\nmax |b | = %8g     max |x0| = %8g", normb , normx0);
	//printf(                "      xsize   = %8g", xsize_);
	printf("\nmax |y0| = %8g     max |z0| = %8g", normy0, normz0);
	//printf(                "      zsize   = %8g", zsize_);

	//---------------------------------------------------------------------
	// Initialize.
	//---------------------------------------------------------------------
	//true   = 1;
	//false  = 0;
	//zn     = zeros(n,1);
	int nb     = n + m;
	int nkkt   = nb;
	int CGitns = 0;
	int inform = 0;
	//---------------------------------------------------------------------  
	//  Only allow scalar d1, d2 for now
	//---------------------------------------------------------------------
	double d1 = D1[0];
	double d2 = D2[0];

	//---------------------------------------------------------------------
	// Grab input options.
	//---------------------------------------------------------------------
	int  maxitn    = params->max_barrier_iterations;
	double featol    = params->feasibility_tolerance;
	double opttol    = params->optimality_tolerance;
	double steptol   = params->step_tolerance;
	int  stepSame  = 1;  /* options.StepSame;   // 1 means stepx == stepz */
	double x0min     = params->x0_min;
	double z0min     = params->z0_min;
	double mu0       = params->mu0;
	int  LSproblem = params->LSproblem;  // See below
	int  LSmethod  = params->LSmethod;   // 1=Cholesky    2=QR    3=LSQR
	RNP::LSQR<double>::Params LSQRParams;
	int  itnlim    = n/2 * std::min(m,n);
	double atol1     = std::numeric_limits<double>::epsilon();  // Initial  atol
	double atol2     = std::numeric_limits<double>::epsilon();  // Smallest atol,unless atol1 is smaller
	double conlim    = 1./std::numeric_limits<double>::epsilon();

	// LSproblem:
	//  1 = dy          2 = dy shifted, DLS
	// 11 = s          12 =  s shifted, DLS    (dx = Ds)
	// 21 = dx
	// 31 = 3x3 system, symmetrized by Z^{1/2}
	// 32 = 2x2 system, symmetrized by X^{1/2}

	//---------------------------------------------------------------------
	// Set other parameters.
	//---------------------------------------------------------------------
	int  kminor    = 0;      // 1 stops after each iteration
	double eta       = 1e-4;   // Linesearch tolerance for "sufficient descent"
	double maxf      = 10;     // Linesearch backtrack limit (function evaluations)
	double maxfail   = 1;      // Linesearch failure limit (consecutive iterations)
	double bigcenter = 1e+3;   // mu is reduced if center < bigcenter.

	// Parameters for LSQR.
	double atolmin   = eps;    // Smallest atol if linesearch back-tracks
	double btol      = 0;      // Should be small (zero is ok)
	double show      = false;  // Controls lsqr iteration log
	/*
	double gamma     = d1->infNorm();
	double delta     = d2->infNorm();
	*/
	double gamma = d1;
	double delta = d2;

	printf("\n\nx0min    = %8g     featol   = %8.1e", x0min, featol);
	printf(                  "      d1max   = %8.1e", gamma);
	printf(  "\nz0min    = %8g     opttol   = %8.1e", z0min, opttol);
	printf(                  "      d2max   = %8.1e", delta);
	printf(  "\nmu0      = %8.1e     steptol  = %8g", mu0  , steptol);
	printf(                  "     bigcenter= %8g"  , bigcenter);

	printf("\n\nLSQR:");
	printf("\natol1    = %8.1e     atol2    = %8.1e", atol1 , atol2 );
	printf(                  "      btol    = %8.1e", btol );
	printf("\nconlim   = %8.1e     itnlim   = %8d"  , conlim, itnlim);
	printf(                  "      show    = %8g"  , show );

	// LSmethod  = 3;  ////// Hardwire LSQR
	// LSproblem = 1;  ////// and LS problem defining "dy".
	/*
	if wait
	printf("\n\nReview parameters... then type "return"\n")
	keyboard
	end
	*/
	if (eta < 0){
		printf("\n\nLinesearch disabled by eta < 0");
	}

	//---------------------------------------------------------------------
	// All parameters have now been set.
	//---------------------------------------------------------------------
	bool useChol = (LSmethod == 1);
	bool useQR   = (LSmethod == 2);
	bool direct  = (LSmethod <= 2 && ifexplicit);

	//---------------------------------------------------------------------
	// Categorize bounds and allow for fixed variables by modifying b.
	//---------------------------------------------------------------------

	int nlow, nupp, nfix;
	size_t *low = new size_t[n];
	size_t *upp = new size_t[n];
	size_t *fix = new size_t[n];
	{ // getBoundTypes(&nlow, &nupp, &nfix, bptrs );
		nlow = 0;
		nupp = 0;
		nfix = 0;
		for(int k = 0; k < n; ++k){
			if(
			bptrs[0][k] = k;
		}
	}

	int nU = n;
	if (nupp ==0) nU = 1;   //Make dummy vectors if no Upper bounds

	//---------------------------------------------------------------------
	//  Get pointers to local copy of model bounds
	//---------------------------------------------------------------------

	double *r1 = new double[m]; RNP::TBLAS::Fill(m, 0., r1, 1);
	double *x1 = new double[n]; RNP::TBLAS::Fill(n, 0., x1, 1);

	if (nfix > 0){
		for (int k=0; k<nfix; k++)
			x1[fix[k]] = bl[fix[k]];
		A(false, m, n, r1, x1, data);
		RNP::TBLAS::Axpy(m, -1., r1, 1, b, 1);
		// At some stage, might want to look at normfix = norm(r1,inf);
	}

	//---------------------------------------------------------------------
	// Scale the input data.
	// The scaled variables are
	//    xbar     = x/beta,
	//    ybar     = y/zeta,
	//    zbar     = z/zeta.
	// Define
	//    theta    = beta*zeta;
	// The scaled function is
	//    phibar   = ( 1   /theta) fbar(beta*xbar),
	//    gradient = (beta /theta) grad,
	//    Hessian  = (beta2/theta) hess.
	//---------------------------------------------------------------------
	double beta = params->xsize;   if (beta==0) beta = 1;   // beta scales b, x.
	double zeta = params->zsize;   if (zeta==0) zeta = 1;   // zeta scales y, z.
	double theta  = beta*zeta;                            // theta scales obj.
	// (theta could be anything, but theta = beta*zeta makes
	// scaled grad = grad/zeta = 1 approximately if zeta is chosen right.)

	for (int k=0; k<nlow; k++)
		bl[low[k]] = bl[low[k]]/beta;
	for (int k=0; k<nupp; k++)
		bu[upp[k]] = bu[upp[k]]/beta;
	d1     = d1 * ( beta/sqrt(theta) );
	d2     = d2 * ( sqrt(theta)/beta );

	double beta2  = beta*beta;
	RNP::TBLAS::Scale(m, 1./beta, b, 1 );
	RNP::TBLAS::Scale(m, 1./zeta, y, 1 );
	RNP::TBLAS::Scale(n, 1./beta, x, 1 );
	RNP::TBLAS::Scale(n, 1./zeta, z, 1 );

	//---------------------------------------------------------------------
	// Initialize vectors that are not fully used if bounds are missing.
	//---------------------------------------------------------------------
	double *rL = new double[n]; RNP::TBLAS::Fill(n, 0., rL, 1);
	double *cL = new double[n]; RNP::TBLAS::Fill(n, 0., cL, 1);
	double *z1 = new double[n]; RNP::TBLAS::Fill(n, 0., z1, 1);
	double *dx1 = new double[n]; RNP::TBLAS::Fill(n, 0., dx1, 1);
	double *dz1 = new double[n]; RNP::TBLAS::Fill(n, 0., dz1, 1);
	double *r2 = new double[n]; RNP::TBLAS::Fill(n, 0., r2, 1);

	double *rU = new double[nU]; RNP::TBLAS::Fill(nU, 0., rU, 1);
	double *cU = new double[nU]; RNP::TBLAS::Fill(nU, 0., cU, 1);
	double *x2 = new double[nU]; RNP::TBLAS::Fill(nU, 0., x2, 1);
	double *z2 = new double[nU]; RNP::TBLAS::Fill(nU, 0., z2, 1);
	double *dx2 = new double[nU]; RNP::TBLAS::Fill(nU, 0., dx2, 1);
	double *dz2 = new double[nU]; RNP::TBLAS::Fill(nU, 0., dz2, 1);

	//---------------------------------------------------------------------
	// Initialize x, y, z, objective, etc.
	//---------------------------------------------------------------------
	double *dx = new double[n]; RNP::TBLAS::Fill(n, 0., dx, 1);
	double *dy = new double[m]; RNP::TBLAS::Fill(m, 0., dy, 1);
	double *Pr = new double[m];
	double *D = new double[n];
	double *w = new double[n];
	double *rhs = new double[m+n];


	//---------------------------------------------------------------------
	// Pull out the element array pointers for efficiency
	//---------------------------------------------------------------------

	for (int k=0; k<nlow; k++){
		x[low[k]]  = std::max( x[low[k]], bl[low[k]]);
		x1[low[k]] = std::max( x[low[k]] - bl[low[k]], x0min  );
		z1[low[k]] = std::max( z[low[k]], z0min  );
	}
	for (int k=0; k<nupp; k++){
		x[upp[k]]  = std::min( x[upp[k]], bu[upp[k]]);
		x2[upp[k]] = std::max(bu[upp[k]] -  x[upp[k]], x0min  );
		z2[upp[k]] = std::max(-z[upp[k]], z0min  );
	}
	//////////////////// Assume hessian is diagonal. //////////////////////

	//  [obj,grad,hess] = feval( Fname, (x*beta) );
	RNP::TBLAS::Scale(n, beta, x, 1 );
	double obj;
	double *grad = new double[n];
	double *H = new double[n];
	phi(n, x, &obj, grad, H, data);
	RNP::TBLAS::Scale(n, 1./beta, x, 1 );

	obj /= theta;                       // Scaled obj.
	for(int k = 0; k < n; ++k)
		grad[k] = grad[k]*(beta/theta) + (d1*d1)*x[k]; // grad includes x regularization.
	for(int k = 0; k < n; ++k)
		H[k]  = H[k]*(beta2/theta) + (d1*d1)      ;     // H    includes x regularization.


	/*---------------------------------------------------------------------
	// Compute primal and dual residuals:
	// r1 =  b - Aprod(x) - d2*d2*y;
	// r2 =  grad - Atprod(y) + z2 - z1;  
	//  rL =  bl - x + x1;
	//  rU =  x + x2 - bu; */
	//---------------------------------------------------------------------
	//  [r1,r2,rL,rU,Pinf,Dinf] = ...
	//      pdxxxresid1( Aname,fix,low,upp, ...
	//                   b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 );
	pdxxxresid1(m, n, A, nlow, nupp, nfix, low, upp, fix, 
		b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2,
		r1, r2, &Pinf, &Dinf, data);
	//---------------------------------------------------------------------
	// Initialize mu and complementarity residuals:
	//    cL   = mu*e - X1*z1.
	//    cU   = mu*e - X2*z2.
	//
	// 25 Jan 2001: Now that b and obj are scaled (and hence x,y,z),
	//              we should be able to use mufirst = mu0 (absolute value).
	//              0.1 worked poorly on StarTest1 with x0min = z0min = 0.1.
	// 29 Jan 2001: We might as well use mu0 = x0min * z0min;
	//              so that most variables are centered after a warm start.
	// 29 Sep 2002: Use mufirst = mu0*(x0min * z0min),
	//              regarding mu0 as a scaling of the initial center.
	//---------------------------------------------------------------------
	//  double mufirst = mu0*(x0min * z0min);
	double mufirst = mu0;   // revert to absolute value
	double mulast  = 0.1 * opttol;
	mulast  = std::min( mulast, mufirst );
	double mu      = mufirst;
	double center,  fmerit;
	pdxxxresid2( mu, nlow, nupp, low,upp, cL, cU, x1, x2,
		z1, z2, &center, &Cinf, &Cinf0 );
		fmerit = pdxxxmerit(m, n, nlow, nupp, low, upp, r1, r2, rL, rU, cL, cU );

	// Initialize other things.

	bool  precon   = true;  
	double PDitns    = 0;
	bool converged = false;
	double atol      = atol1;
	atol2     = std::max( atol2, atolmin );
	atolmin   = atol2;
	//  pdDDD2    = d2;    // Global vector for diagonal matrix D2

	//  Iteration log.

	double stepx   = 0;
	double stepz   = 0;
	int nf      = 0;
	int itncg   = 0;
	int nfail   = 0;

	printf("\n\nItn   mu   stepx   stepz  Pinf  Dinf");
	printf("  Cinf   Objective    nf  center");
	if (direct) {
		printf("\n");
	}else{ 
		printf("  atol   solver   Inexact\n");
	}

	double regx = d1*RNP::TBLAS::Norm2(n, x, 1);
	double regy = d2*RNP::TBLAS::Norm2(m, y, 1);
	//  regterm = twoNorm(d1.*x)^2  +  norm(d2.*y)^2;
	double regterm = regx*regx + regy*regy;
	double objreg  = obj  +  0.5*regterm;
	double objtrue = objreg * theta;

	printf("\n%3g                     ", PDitns        );
	printf("%6.1f%6.1f" , log10(Pinf ), log10(Dinf));
	printf("%6.1f%15.7e", log10(Cinf0), objtrue    );
	printf("   %8.1f\n"   , center                   );

	//---------------------------------------------------------------------
	// Main loop.
	//---------------------------------------------------------------------
	Info info;
	Outfo outfo;
	while(PDitns < maxitn) {
		++PDitns;

		// 31 Jan 2001: Set atol according to progress, a la Inexact Newton.
		// 07 Feb 2001: 0.1 not small enough for Satellite problem.  Try 0.01.
		// 25 Apr 2001: 0.01 seems wasteful for Star problem.
		//              Now that starting conditions are better, go back to 0.1.

		double r3norm = std::max(Pinf,   std::max(Dinf,  Cinf));
		atol   = std::min(atol,  r3norm*0.1);
		atol   = std::max(atol,  atolmin   );
		info.r3norm = r3norm;

		//-------------------------------------------------------------------
		//  Define a damped Newton iteration for solving f = 0,
		//  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
		//  to obtain the system
		//
		//     [-H2  A"  ] [ dx ] = [ w ],   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
		//     [ A   D2^2] [ dy ] = [ r1]    w  = r2 - X1inv(cL + Z1 rL)
		//                                           + X2inv(cU + Z2 rU),
		//
		//  which is equivalent to the least-squares problem
		//
		//     min || [ D A"]dy  -  [  D w   ] ||,   D = H2^{-1/2}.         (*)
		//         || [  D2 ]       [D2inv r1] ||
		//-------------------------------------------------------------------
		for (int k=0; k<nlow; k++)
			H[low[k]]  = H[low[k]] + z1[low[k]]/x1[low[k]];
		for (int k=0; k<nupp; k++)
			H[upp[k]]  = H[upp[k]] + z2[upp[k]]/x2[upp[k]];
		w = r2;
		for (int k=0; k<nlow; k++)
			w[low[k]]  = w[low[k]] - (cL[low[k]] + z1[low[k]]*rL[low[k]])/x1[low[k]];
		for (int k=0; k<nupp; k++)
			w[upp[k]]  = w[upp[k]] + (cU[upp[k]] + z2[upp[k]]*rU[upp[k]])/x2[upp[k]];

		if (LSproblem == 1){
			//-----------------------------------------------------------------
			//  Solve (*) for dy.
			//-----------------------------------------------------------------
			for(int k = 0; k < n; ++k)
				H[k]      = 1.0/H[k];    // H is now Hinv (NOTE!)
			for (int k=0; k<nfix; k++)  
				H[fix[k]] = 0;
			for (int k=0; k<n; k++)
				D[k]= sqrt(H[k]);
			//thisLsqr.borrowDiag1(D);
			//thisLsqr.diag2_ = d2;

			if (direct){
				// Omit direct option for now
			}else {// Iterative solve using LSQR.
				//rhs     = [ D.*w; r1./d2 ];
				for (int k=0; k<n; k++)
				rhs[k] = D[k]*w[k];
				for (int k=0; k<m; k++)
				rhs[n+k] = r1[k]*(1.0/d2);
				double damp    = 0;

				//if (precon){    // Construct diagonal preconditioner for LSQR
				//	matPrecon(d2, Pr, D);
				//}  

				//  New version of lsqr

				int istop, itn;
				RNP::TBLAS::Fill(m, 0., dy, 1);
				show = false;
				info.atolmin = atolmin;
				info.r3norm  = fmerit;  // Must be the 2-norm here.

				do_lsqr(m, n, A, D, d2, rhs, damp, atol, btol, conlim, itnlim, show, info, dy , &istop, &itncg, &outfo, precon, Pr, data);
				/*
				if (precon)
					dy = dy*Pr;

				if (!precon && itncg>999999)
					precon=true;
				*/
				if (istop == 3  ||  istop == 7 )  // conlim or itnlim
					printf("\n    LSQR stopped early:  istop = //%d", istop);


				atolold   = outfo.atolold;
				atol      = outfo.atolnew;
				r3ratio   = outfo.r3ratio;
			}// LSproblem 1

			//      grad      = pdxxxmat( Aname,2,m,n,dy );   // grad = A"dy
			RNP::TBLAS::Fill(n, 0., grad, 1);
			A(true, m, n, grad, dy, data);
			for (int k=0; k<nfix; k++)
				grad[fix[k]] = 0;                            // grad is a work vector
			for(int k = 0; k < n; ++k){
				dx[k] = H[k] * (grad[k] - w[k]);
			}
		}else{
			perror( "This LSproblem not yet implemented\n" );
		}
		//-------------------------------------------------------------------

		CGitns += itncg;

		//-------------------------------------------------------------------
		// dx and dy are now known.  Get dx1, dx2, dz1, dz2.
		//-------------------------------------------------------------------
		for (int k=0; k<nlow; k++){    
			dx1[low[k]] = - rL[low[k]] + dx[low[k]];
			dz1[low[k]] =  (cL[low[k]] - z1[low[k]]*dx1[low[k]]) / x1[low[k]];
		}
		for (int k=0; k<nupp; k++){    
			dx2[upp[k]] = - rU[upp[k]] - dx[upp[k]];
			dz2[upp[k]] =  (cU[upp[k]] - z2[upp[k]]*dx2[upp[k]]) / x2[upp[k]];
		}
		//-------------------------------------------------------------------
		// Find the maximum step.
		//--------------------------------------------------------------------
		double stepx1 = pdxxxstep(n, nlow, low, x1, dx1 );
		double stepx2 = inf;
		if (nupp > 0)
			stepx2 = pdxxxstep(n, nupp, upp, x2, dx2 );
		double stepz1 = pdxxxstep(n,  z1     , dz1      );
		double stepz2 = inf;
		if (nupp > 0)
			stepz2 = pdxxxstep(n, z2     , dz2      );
		double stepx  = std::min( stepx1, stepx2 );
		double stepz  = std::min( stepz1, stepz2 );
		stepx  = std::min( steptol*stepx, 1.0 );
		stepz  = std::min( steptol*stepz, 1.0 );
		if (stepSame){                   // For NLPs, force same step
			stepx = std::min( stepx, stepz );   // (true Newton method)
			stepz = stepx;
		}

		//-------------------------------------------------------------------
		// Backtracking linesearch.
		//-------------------------------------------------------------------
		bool fail     =  true;
		nf       =  0;

		while (nf < maxf){
			nf      = nf + 1;
			RNP::TBLAS::Axpy(n, stepx, dx, 1, x, 1);
			RNP::TBLAS::Axpy(m, stepz, dy, 1, y, 1);
			for (int k=0; k<nlow; k++){    
				x1[low[k]] = x1[low[k]]  +  stepx * dx1[low[k]];
				z1[low[k]] = z1[low[k]]  +  stepz * dz1[low[k]];
			}
			for (int k=0; k<nupp; k++){    
				x2[upp[k]] = x2[upp[k]]  +  stepx * dx2[upp[k]];
				z2[upp[k]] = z2[upp[k]]  +  stepz * dz2[upp[k]];
			}
			//      [obj,grad,hess] = feval( Fname, (x*beta) );
			RNP::TBLAS::Scale(n, beta, x, 1 );
			phi(n, x, &obj, grad, H, data);
			RNP::TBLAS::Scale(n, 1./beta, x, 1 );

			obj        /= theta;
			for(int k = 0; k < n; ++k)
				grad[k]       = grad[k]*(beta /theta)  +  d1*d1*x[k];
			for(int k = 0; k < n; ++k)
				H[k]          = H[k]*(beta2/theta)  +  d1*d1;

			//      [r1,r2,rL,rU,Pinf,Dinf] = ...
			pdxxxresid1(m, n, A, nlow, nupp, nfix, low, upp, fix,
				b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,
				y,z1,z2, r1, r2, &Pinf, &Dinf, data );
			//double center, Cinf, Cinf0;
			//      [cL,cU,center,Cinf,Cinf0] = ...
			pdxxxresid2( mu, nlow, nupp, low,upp,cL,cU,x1,x2,z1,z2,
				&center, &Cinf, &Cinf0);
			double fmeritnew = pdxxxmerit(m, n, nlow, nupp, low,upp,r1,r2,rL,rU,cL,cU );
			double step      = std::min( stepx, stepz );

			if (fmeritnew <= (1 - eta*step)*fmerit){
				fail = false;
				break;
			}

			// Merit function didn"t decrease.
			// Restore variables to previous values.
			// (This introduces a little error, but save lots of space.)

			RNP::TBLAS::Axpy(n, -stepx, dx, 1, x, 1);
			RNP::TBLAS::Axpy(m, -stepz, dy, 1, y, 1);
			for (int k=0; k<nlow; k++){    
				x1[low[k]] = x1[low[k]]  -  stepx * dx1[low[k]];
				z1[low[k]] = z1[low[k]]  -  stepz * dz1[low[k]];
			}
			for (int k=0; k<nupp; k++){    
				x2[upp[k]] = x2[upp[k]]  -  stepx * dx2[upp[k]];
				z2[upp[k]] = z2[upp[k]]  -  stepz * dz2[upp[k]];
			}
			// Back-track.
			// If it"s the first time,
			// make stepx and stepz the same.

			if (nf == 1 && stepx != stepz){
				stepx = step;
			}else if (nf < maxf){
				stepx = stepx/2;
			}
			stepz = stepx;
		}

		if (fail){
			printf("\n     Linesearch failed (nf too big)");
			nfail += 1;
		}else{
			nfail = 0;
		}

		//-------------------------------------------------------------------
		// Set convergence measures.
		//--------------------------------------------------------------------
		regx = d1*RNP::TBLAS::Norm2(n, x, 1);
		regy = d2*RNP::TBLAS::Norm2(m, y, 1);
		regterm = regx*regx + regy*regy;
		objreg  = obj  +  0.5*regterm;
		objtrue = objreg * theta;

		bool primalfeas    = Pinf  <=  featol;
		bool dualfeas      = Dinf  <=  featol;
		bool complementary = Cinf0 <=  opttol;
		bool enough        = PDitns>=       4;  // Prevent premature termination.
		bool converged     = primalfeas  &  dualfeas  &  complementary  &  enough;

		//-------------------------------------------------------------------
		// Iteration log.
		//-------------------------------------------------------------------
		char str1[100],str2[100],str3[100],str4[100],str5[100];
		sprintf(str1, "\n%3g%5.1f" , PDitns      , log10(mu)   );
		sprintf(str2, "%8.5f%8.5f" , stepx       , stepz       );
		if (stepx < 0.0001 || stepz < 0.0001){
			sprintf(str2, " %6.1e %6.1e" , stepx       , stepz       );
		}

		sprintf(str3, "%6.1f%6.1f" , log10(Pinf) , log10(Dinf));
		sprintf(str4, "%6.1f%15.7e", log10(Cinf0), objtrue     );
		sprintf(str5, "%3d%8.1f"   , nf          , center      );
		if (center > 99999){
			sprintf(str5, "%3d%8.1e"   , nf          , center      );
		}
		printf("%s%s%s%s%s", str1, str2, str3, str4, str5);
		if (direct){
			// relax
		}else{
			printf(" %5.1f%7d%7.3f", log10(atolold), itncg, r3ratio);
		}
		//-------------------------------------------------------------------
		// Test for termination.
		//-------------------------------------------------------------------
		if (kminor){
			printf( "\nStart of next minor itn...\n");
			//      keyboard;
		}

		if (converged){
			printf("\n   Converged");
			break;
		}else if (PDitns >= maxitn){
			printf("\n   Too many iterations");
			inform = 1;
			break;
		}else if (nfail  >= maxfail){
			printf("\n   Too many linesearch failures");
			inform = 2;
			break;
		}else{

			// Reduce mu, and reset certain residuals.

			double stepmu  = std::min( stepx , stepz   );
			stepmu  = std::min( stepmu, steptol );
			double muold   = mu;
			mu      = mu   -  stepmu * mu;
			if (center >= bigcenter)
				mu = muold;  

			// mutrad = mu0*(sum(Xz)/n); // 24 May 1998: Traditional value, but
			// mu     = CoinMin(mu,mutrad ); // it seemed to decrease mu too much.

			mu      = std::max(mu,mulast);  // 13 Jun 1998: No need for smaller mu.
			//      [cL,cU,center,Cinf,Cinf0] = ...
			pdxxxresid2( mu, nlow, nupp, low, upp, cL, cU, x1, x2, z1, z2,
				&center, &Cinf, &Cinf0 );
				fmerit = pdxxxmerit(m, n, nlow, nupp, low,upp,r1,r2,rL,rU,cL,cU );

			// Reduce atol for LSQR (and SYMMLQ).
			// NOW DONE AT TOP OF LOOP.

			atolold = atol;
			// if atol > atol2
			//   atolfac = (mu/mufirst)^0.25;
			//   atol    = CoinMax( atol*atolfac, atol2 );
			// end

			// atol = CoinMin( atol, mu );     // 22 Jan 2001: a la Inexact Newton.
			// atol = CoinMin( atol, 0.5*mu ); // 30 Jan 2001: A bit tighter

			// If the linesearch took more than one function (nf > 1),
			// we assume the search direction needed more accuracy
			// (though this may be true only for LPs).
			// 12 Jun 1998: Ask for more accuracy if nf > 2.
			// 24 Nov 2000: Also if the steps are small.
			// 30 Jan 2001: Small steps might be ok with warm start.
			// 06 Feb 2001: Not necessarily.  Reinstated tests in next line.

			if (nf > 2  ||  std::min( stepx, stepz ) <= 0.01)
				atol = atolold*0.1;
		}
		//---------------------------------------------------------------------
		// End of main loop.
		//---------------------------------------------------------------------
	}


	for (int k=0; k<nfix; k++) x[fix[k]] = bl[fix[k]];
	RNP::TBLAS::Copy(n, z1, 1, z, 1);
	if (nupp > 0){
		RNP::TBLAS::Axpy(n, -1., z2, 1, z, 1);
	}
	printf("\n\nmax |x| =%10.3f", x[RNP::TBLAS::MaximumIndex(n, x, 1)] );
	printf("    max |y| =%10.3f", y[RNP::TBLAS::MaximumIndex(m, y, 1)] );
	printf("    max |z| =%10.3f", z[RNP::TBLAS::MaximumIndex(n, z, 1)] );
	printf("   scaled");

	// Unscale x, y, z.
	RNP::TBLAS::Scale(n, beta, x, 1);
	RNP::TBLAS::Scale(m, beta, y, 1);
	RNP::TBLAS::Scale(n, beta, z, 1); 

	printf(  "\nmax |x| =%10.3f", x[RNP::TBLAS::MaximumIndex(n, x, 1)] );
	printf("    max |y| =%10.3f", y[RNP::TBLAS::MaximumIndex(m, y, 1)] );
	printf("    max |z| =%10.3f", z[RNP::TBLAS::MaximumIndex(n, z, 1)] );
	printf(" unscaled\n");

	char str1[100],str2[100];
	sprintf(str1, "\nPDitns  =%10g", PDitns );
	sprintf(str2, "itns =%10d", CGitns );
	//  printf( [str1 " " solver str2] );

	//-----------------------------------------------------------------------
	// End function pdco.m
	//-----------------------------------------------------------------------

	// Print distribution
	float thresh[9]={ 0.00000001, 0.0000001, 0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1.00001};
	int counts[9]={0};
	for (int ij=0; ij<n; ij++){
		for (int j=0; j<9; j++){
			if(x[ij] < thresh[j]){
				counts[j] += 1;
				break;
			}
		}
	}
	printf ("Distribution of Solution Values\n");
	for (int j=8; j>1; j--)
		printf(" %f  to  %f %d\n",thresh[j-1],thresh[j],counts[j]);
	printf("   Less than   %f %d\n",thresh[2],counts[0]);

	delete [] bptrs[0];
	delete [] bptrs[1];
	delete [] bptrs[2];
	delete [] r1;
	delete [] x1;
	delete [] rL;
	delete [] cL;
	delete [] z1;
	delete [] dx1;
	delete [] dz1;
	delete [] r2;

	delete [] rU;
	delete [] cU;
	delete [] x2;
	delete [] z2;
	delete [] dx2;
	delete [] dz2;

	delete [] dx;
	delete [] dy;
	delete [] Pr;
	delete [] D;
	delete [] w;
	delete [] rhs;
	
	delete [] grad;
	delete [] H;
	return 0;
}


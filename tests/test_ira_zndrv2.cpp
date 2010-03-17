#include <cstdio>
#include <cstring>
#include <complex>
#include "TBLAS.h"
#include "IRA.h"

// Expected results:
//
// Ritz values (Real, Imag) and relative residuals
// -----------------------------------------------
//               Col   1       Col   2       Col   3
//  Row   1:    3.48720E+01   2.04672E-14   2.24418E-13
//  Row   2:    6.44326E+01  -3.79074E-14   1.35536E-13
//  Row   3:    1.13669E+02   2.42490E-13   7.86603E-14
//  Row   4:    1.82532E+02  -7.55066E-13   5.54891E-14

static const size_t maxn = 256;
static const size_t maxncv = 30;
static const size_t ldv = maxn;
static const double rho = 10.;

struct OpData{
	std::complex<double> dd[maxn], dl[maxn], du[maxn], du2[maxn];
	size_t ipiv[maxn];
};


static inline double _abs1(const std::complex<double> &z){
	return fabs(z.real()) + fabs(z.imag());
}

void av_(size_t n, const std::complex<double> *v, std::complex<double> *w){
	const double h = 1./(n+1);
	const double h2 = h*h;
	const double s = rho/2.0;
	const double dd = 2.0/h2;
	const double dl = -1./h2 - s/h;
	const double du = -1./h2 + s/h;
	w[0] = dd*v[0] + du*v[1];
	for(size_t j = 1; j < n-1; ++j){
		w[j] = dl*v[j-1] + dd*v[j] + du*v[j+1];
	}
	w[n-1] = dl*v[n-2] + dd*v[n-1];
}
void op_(size_t n, const std::complex<double> &sigma, const std::complex<double> *v, std::complex<double> *w, void *data){
	OpData *op_data = (OpData*)data;
	RNP::TBLAS::Copy(n, v, 1, w, 1);
	// zgttrs('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr)
	{
		// Solve L*x = b.
		for(size_t i = 0; i < n-1; ++i){
			if(op_data->ipiv[i] == i){
				w[i+1] -= op_data->dl[i]*w[i];
			}else{
				std::complex<double> temp = w[i];
				w[i] = w[i+1];
				w[i+1] = temp - op_data->dl[i]*w[i];
			}
		}
		// Solve U*x = b.
		w[n-1] /= op_data->dd[n-1];
		if(n > 1){
			w[n-2] = (w[n-2]-op_data->du[n-2]*w[n-1]) / op_data->dd[n-2];
		}
		{size_t i = n-2;
		while(i --> 0){
			w[i] = (w[i]-op_data->du[i]*w[i+1] - op_data->du2[i]*w[i+2]) / op_data->dd[i];
		}}
	}
}

int main(int argc, char *argv[]){
// 
//      Simple program to illustrate the idea of reverse communication
//      in shift-invert mode for a standard complex nonsymmetric eigenvalue 
//      problem.
// 
//      We implement example two of ex-complex.doc in DOCUMENTS directory
// 
// \Example-2
//      ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
//          where A is derived from the central difference discretization
//          of the 1-dimensional convection-diffusion operator
//                    (d^2u/dx^2) + rho*(du/dx)
//          on the interval [0,1] with zero Dirichlet boundary condition.
//      ... The shift sigma is a complex number.
// 
//      ... OP = inv[A-sigma*I] and  B = I.
// 
//      ... Use mode 3 of ZNAUPD.
// 
// \BeginLib
// 
// \Routines called:
//      znaupd  ARPACK reverse communication interface routine.
//      zneupd  ARPACK routine that returns Ritz values and (optionally)
//              Ritz vectors.
//      zgttrf  LAPACK tridiagonal factorization routine.
//      zgttrs  LAPACK tridiagonal solve routine.
//      dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
//      zaxpy   Level 1 BLAS that computes y <- alpha*x+y.
//      zcopy   Level 1 BLAS that copies one vector to another.
//      dznrm2  Level 1 BLAS that computes the norm of a vector.
//      av      Matrix vector multiplication routine that computes A*x.
// 
// \Author
//      Richard Lehoucq
//      Danny Sorensen
//      Chao Yang
//      Dept. of Computational &
//      Applied Mathematics
//      Rice University
//      Houston, Texas
// 
// \SCCS Information: @(#)
//  FILE: ndrv2.F   SID: 2.2   DATE OF SID: 4/22/96   RELEASE: 2
// 
// \Remarks
//      1. None
// 
// \EndLib
// --------------------------------------------------------------------------
// 
	std::complex<double> ax[maxn], d[maxncv], v[ldv*maxncv];
	OpData op_data;
	double rd[maxncv*3];

	size_t nx = 10;
	size_t n = nx*nx;
	size_t nev = 4;
	size_t ncv = 20;
	
	std::complex<double> sigma(0.0);
	
	// Construct C = A - SIGMA*I, factor C in complex
	// arithmetic (using LAPACK subroutine zgttrf). The
	// matrix A is chosen to be the tridiagonal matrix
	// derived from standard central difference of the
	// 1-d convection diffusion operator - u" + rho*u' on
	// the interval [0, 1] with zero Dirichlet boundary
	// condition.
	const double h = 1.0/(n+1.0);
	const double h2 = h*h;
	const double s = 0.5*rho;
	const double s1 = -1.0/h2 - s/h;
	const std::complex<double> s2 = 2.0/h2 - sigma;
	const double s3 = -1.0/h2 + s/h;
	for(size_t j = 0; j < n-1; ++j){
		op_data.dl[j] = s1;
		op_data.dd[j] = s2;
		op_data.du[j] = s3;
	}
	op_data.dd[n-1] = s2;
	// zgttrf(n, dl, dd, du, du2, ipiv, ierr)
	{
		for(size_t i = 0; i < n; ++i){
			op_data.ipiv[i] = i;
		}
		for(size_t i = 0; i < n-2; ++i){
			op_data.du2[i] = 0;
		}
		for(size_t i = 0; i < n-2; ++i){
			if(_abs1(op_data.dd[i]) >= _abs1(op_data.dl[i])){
				// No row interchange required, eliminate dl[i]
				if(_abs1(op_data.dd[i]) != 0){
					std::complex<double> fact = op_data.dl[i] / op_data.dd[i];
					op_data.dl[i] = fact;
					op_data.dd[i+1] -= fact * op_data.du[i];
				}
			}else{
				// Interchange rows i and i+1, eliminate dl[i]
				std::complex<double> fact = op_data.dd[i] / op_data.dl[i];
				op_data.dd[i] = op_data.dl[i];
				op_data.dl[i] = fact;
				std::complex<double> temp = op_data.du[i];
				op_data.du[i] = op_data.dd[i+1];
				op_data.dd[i+1] = temp - fact*op_data.dd[i+1];
				op_data.du2[i] = op_data.du[i+1];
				op_data.du[i+1] *= -fact;
				op_data.ipiv[i] = i+1;
			}
		}
		if(n > 1){
			size_t i = n - 2;
			if(_abs1(op_data.dd[i]) >= _abs1(op_data.dl[i])){
				if(_abs1(op_data.dd[i]) != 0){
					std::complex<double> fact = op_data.dl[i] / op_data.dd[i];
					op_data.dl[i] = fact;
					op_data.dd[i+1] -= fact*op_data.du[i];
				}
			}else{
				std::complex<double> fact = op_data.dd[i] / op_data.dl[i];
				op_data.dd[i] = op_data.dl[i];
				op_data.dl[i] = fact;
				std::complex<double> temp = op_data.du[i];
				op_data.du[i] = op_data.dd[i+1];
				op_data.dd[i+1] = temp - fact*op_data.dd[i+1];
				op_data.ipiv[i] = i + 1;
			}
		}
	}
	
	int nconv = RNP::IRA::ShiftInvert(
		n, sigma, &op_, NULL,
		nev, ncv, &RNP::IRA::LargestMagnitude,
		d, v, ldv,
		NULL,
		NULL,
		(void*)&op_data);
	if(nconv > 0){
		for(int j = 0; j < nconv; ++j) {
			// Compute the residual norm
			//
			//   ||  A*x - lambda*x ||
			//
			// for the NCONV accurately computed eigenvalues and eigenvectors.
			// (iparam(5) indicates how many are accurate to the requested tolerance)

			av_(n, &v[0+j*ldv], ax);
			RNP::TBLAS::Axpy(n, -d[j], &v[0+j*ldv], 1, ax, 1);
			rd[j+0*maxncv] = d[j].real();
			rd[j+1*maxncv] = d[j].imag();
			rd[j+2*maxncv] = RNP::TBLAS::Norm2(n, ax, 1);
			rd[j+2*maxncv] /= std::abs(d[j]);
		}
		
		// Display computed residuals.
		printf("Residuals:\n");
		for(int i = 0; i < nconv; ++i){
			for(size_t j = 0; j < 3; ++j){
				printf("\t%e", rd[i+j*maxncv]);
			}printf("\n");
		}
	}else{
		printf("IRA::ShiftInvert returned %d\n", nconv);
	}
	return 0;
}


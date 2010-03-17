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
//  Row   1:    7.16197E+02   1.02958E+03   9.01377E-15
//  Row   2:    7.16197E+02  -1.02958E+03   9.38682E-15
//  Row   3:    6.87583E+02   1.02958E+03   2.50909E-15
//  Row   4:    6.87583E+02  -1.02958E+03   2.82193E-15


void tv_(size_t nx, const std::complex<double> *x, std::complex<double> *y){
	// Compute the matrix vector multiplication y<---T*x
	// where T is a nx by nx tridiagonal matrix with DD on the
	// diagonal, DL on the subdiagonal, and DU on the superdiagonal

    double h = 1./(nx+1);
    double h2 = h*h;
    double dd = 4./h2;
	double dl  = -1./h2 - 0.5 *100./h;
	double du  = -1./h2 + 0.5 *100./h;

	y[0] = dd*x[0] + du*x[1];
	for(size_t j = 1; j < nx-1; ++j){
		y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
	}
	y[nx-1] =  dl*x[nx-2] + dd*x[nx-1];
}

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the convection-diffusion operator */
/*     discretized using centered difference. */

void av_(size_t n, const std::complex<double> *v, std::complex<double> *w, void *data){
	size_t nx = *(size_t*)data;
    size_t lo;
    extern void tv_(size_t, const std::complex<double> *, std::complex<double> *);

	// Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
	// tridiagonal matrix

	//              | T -I          |
	//              |-I  T -I       |
	//         OP = |   -I  T       |
	//              |        ...  -I|
	//              |           -I T|

	// derived from the standard central difference  discretization
	// of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
	// with zero boundary condition. */

	// The subroutine TV is called to computed y<---T*x.

    double h2 = 1./((nx + 1) * (nx + 1));

    tv_(nx, v, w);
    RNP::TBLAS::Axpy(nx, -1./h2, &v[nx], 1, w, 1);

    for (size_t j = 1; j < nx-1; ++j) {
		lo = j * nx;
		tv_(nx, &v[lo], &w[lo]);
		RNP::TBLAS::Axpy(nx, -1./h2, &v[lo - nx], 1, &w[lo], 1);
		RNP::TBLAS::Axpy(nx, -1./h2, &v[lo + nx], 1, &w[lo], 1);
    }

    lo = (nx - 1) * nx;
    tv_(nx, &v[lo], &w[lo]);
	RNP::TBLAS::Axpy(nx, -1./h2, &v[lo - nx], 1, &w[lo], 1);
}

int main(int argc, char *argv[]){
	static const size_t maxn = 256;
	static const size_t maxncv = 30;
	static const size_t ldv = maxn;
	
	std::complex<double> ax[maxn], d[maxncv], v[ldv*maxncv];
	double rd[maxncv*3];

/*     Example program to illustrate the idea of reverse communication */
/*     for a standard complex nonsymmetric eigenvalue problem. */

/*     We implement example one of ex-complex.doc in DOCUMENTS directory */

/* \Example-1 */
/*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
/*         where A is obtained from the standard central difference */
/*         discretization of the convection-diffusion operator */
/*                 (Laplacian u) + rho*(du / dx) */
/*         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary */
/*         condition. */

/*     ... OP = A  and  B = I. */

/*     ... Assume "call av (nx,x,y)" computes y = A*x */

/*     ... Use mode 1 of ZNAUPD. */

/* \BeginLib */

/* \Routines called */
/*     znaupd  ARPACK reverse communication interface routine. */
/*     zneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     av      Matrix vector multiplication routine that computes A*x. */
/*     tv      Matrix vector multiplication routine that computes T*x, */
/*             where T is a tridiagonal matrix.  It is used in routine */
/*             av. */

/* \Author */
/*     Richard Lehoucq */
/*     Danny Sorensen */
/*     Chao Yang */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: ndrv1.F   SID: 2.2   DATE OF SID: 4/22/96   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */
/* --------------------------------------------------------------------------- */

	// Define maximum dimensions
	// for all arrays.
	// MAXN:   Maximum dimension
	//         of the A allowed.
	// MAXNEV: Maximum NEV allowed
	// MAXNCV: Maximum NCV allowed

	// The number NX is the number of interior points
	// in the discretization of the 2-dimensional
	// convection-diffusion operator on the unit
	// square with zero Dirichlet boundary condition.
	// The number N(=NX*NX) is the dimension of the
	// matrix.  A standard eigenvalue problem is
	// solved (BMAT = 'I').  NEV is the number of
	// eigenvalues to be approximated.  The user can
	// modify NX, NEV, NCV, WHICH to solve problems of
	// different sizes, and to get different parts of
	// the spectrum.  However, The following
	// conditions must be satisfied:
	//                   N <= MAXN
	//                 NEV <= MAXNEV
	//           NEV + 2 <= NCV <= MAXNCV

	size_t nx = 10;
	size_t n = nx*nx;
	size_t nev = 4;
	size_t ncv = 20;
	
	int nconv = RNP::IRA::Regular(
		n, &av_, NULL,
		nev, ncv, &RNP::IRA::LargestMagnitude,
		d, v, ldv,
		NULL,
		NULL,
		(void*)&nx);
	if(nconv > 0){
		for(int j = 0; j < nconv; ++j) {
			// Compute the residual norm
			//
			//   ||  A*x - lambda*x ||
			//
			// for the NCONV accurately computed eigenvalues and eigenvectors.
			// (iparam(5) indicates how many are accurate to the requested tolerance)

			av_(nx, &v[0+j*ldv], ax, (void*)&nx);
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
		printf("IRA::Regular returned %d\n", nconv);
	}
	return 0;
}


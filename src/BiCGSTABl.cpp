#include <cstddef>
#include <cmath>
#include <complex>
#include "TBLAS.h"
#include "TLASupport.h"
#include "IterativeLinearSolvers.h"

// According to the author, Gerard Sleijpen:
// The code is distributed under the terms of the GNU General Public
// License (version 2 of the License, or any later version) as
// published by the Free Software Foundation.
//
// (However, much of the code has been rewritten.)

void RNP::BiCGSTABl_alloc(size_t n, const BiCGSTABl_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	const size_t l = params.l;
	workspace->zwork = new std::complex<double>[n*(3+2*(l+1)) + (l+1)*(3+2*(l+1))]; // the first term is work, second is rwork, if this were a purely real formulation
	workspace->iwork = new size_t[l+1];
}
void RNP::BiCGSTABl_free(size_t n, const BiCGSTABl_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	if(NULL != workspace->zwork){ delete [] workspace->zwork; workspace->zwork = NULL; }
	if(NULL != workspace->iwork){ delete [] workspace->iwork; workspace->iwork = NULL; }
}

int RNP::BiCGSTABl(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	RNP::LinearSolverParams &ls_params,
	const RNP::BiCGSTABl_params &params,
	// Optional parameters
	RNP::LinearSolverWorkspace *workspace,
	//std::complex<double> *work, // if not null, should be length (3*n + 2*s)*(s+1)
	//size_t *iwork, // if not null, should be length params.s
	void *Adata,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*),
	void *Pdata
){
	const size_t l = params.l;

	RNP::LinearSolverWorkspace *_work = workspace;
	if(NULL == workspace){
		_work = new RNP::LinearSolverWorkspace();
		RNP::BiCGSTABl_alloc(n, params, _work);
	}
	size_t *iwork = _work->iwork;
	std::complex<double> *rwork = _work->zwork + n*(3+2*(l+1));
	int info = 0;

	std::complex<double> *rr = _work->zwork;
	std::complex<double> *r = rr + n;
	std::complex<double> *u = r + n*(l+1);
	std::complex<double> *xp = u + n*(l+1);
	std::complex<double> *bp = xp + n;

	std::complex<double> *z = rwork;
	std::complex<double> *zz = z + (l+1)*(l+1);
	std::complex<double> *y0 = zz + (l+1)*(l+1);
	std::complex<double> *yl = y0 + (l+1);
	std::complex<double> *y = yl + (l+1);
	
	// Set initial x
	if(!ls_params.x_initialized){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
	}
	// Initialize first residual
	A(x, r, Adata);
	for(size_t i = 0; i < n; ++i){
		r[i] = b[i] - r[i];
	}
	if(NULL != P){
		P(r, Pdata);
	}

	size_t nmv = 0;

	RNP::TBLAS::Copy(n, r, 1, rr, 1);
	RNP::TBLAS::Copy(n, r, 1, bp, 1);
	RNP::TBLAS::Copy(n, x, 1, xp, 1);
	RNP::TBLAS::SetMatrix<'N'>(n, 1, 0., 0., x, 1);
	const double rnrm0 = RNP::TBLAS::Norm2(n, r, 1);
	double rnrm = rnrm0;

	double mxnrmx = rnrm0;
	double mxnrmr = rnrm0;
	
	std::complex<double> alpha = 0.;
	std::complex<double> omega = 1.;
	std::complex<double> sigma = 1.;
	std::complex<double> rho0 = 1.;
	size_t iter = 0;
	if(0 == ls_params.max_mult){ ls_params.max_mult = (size_t)-1; }
	while(rnrm > ls_params.tol*rnrm0  &&  nmv < ls_params.max_mult){
		++iter;
		// The BiCG part
		rho0 = -omega*rho0;
		for(size_t k = 0; k < l; ++k){
			std::complex<double> rho1 = RNP::TBLAS::ConjugateDot(n, rr, 1, &r[0+k*n], 1);
			if (rho0 == 0.) {
				return 2;
			}
			std::complex<double> beta = alpha*(rho1/rho0);
			rho0 = rho1;
			for(size_t j = 0; j <= k; ++j){
				for(size_t i = 0; i < n; ++i){
					u[i+j*n] = r[i+j*n] - beta*u[i+j*n];
				}
			}
			A(&u[0+k*n], &u[0+(k+1)*n], Adata);
			if(NULL != P){
				P(&u[0+(k+1)*n], Pdata);
			}
			++nmv;
			sigma = RNP::TBLAS::ConjugateDot(n, rr, 1, &u[0+(k+1)*n], 1);
			if(sigma == 0.){
				return 2;
			}
			alpha = rho1/sigma;
			RNP::TBLAS::Axpy(n, alpha, &u[0+0*n], 1, x, 1);
			for(size_t j = 0; j <= k; ++j){
				RNP::TBLAS::Axpy(n, -alpha, &u[0+(j+1)*n], 1, &r[0+j*n], 1);
			}
			A(&r[0+k*n], &r[0+(k+1)*n], Adata);
			if(NULL != P){
				P(&r[0+(k+1)*n], Pdata);
			}
			++nmv;
			rnrm = RNP::TBLAS::Norm2(n, r, 1);
			if(rnrm > mxnrmx){ mxnrmx = rnrm; }
			if(rnrm > mxnrmr){ mxnrmr = rnrm; }
		}
		// The convex polynomial part
		// Z = R'R
		for(size_t i = 0; i <= l; ++i){
			RNP::TBLAS::MultMV<'C'>(n, l+1-i, 1., &r[0+i*n], n, &r[0+i*n], 1, 0., &z[i+i*(l+1)], 1);
			if (l-i != 0) {
				RNP::TBLAS::Copy(l-i, &z[i+1+i*(l+1)], 1, &z[i+(i+1)*(l+1)], l+1);
				RNP::TBLAS::Conjugate(l-i, &z[i+(i+1)*(l+1)], l+1);
			}
		}
		RNP::TBLAS::CopyMatrix<'A'>(l+1, l+1, z, l+1, zz, l+1);
		RNP::TLASupport::LUDecomposition(l-1, l-1, &zz[1+1*(l+1)], l+1, iwork);
		// tilde r0 and tilde rl (small vectors)
		y0[0] = -1.;
		RNP::TBLAS::Copy(l-1, &z[1+0*(l+1)], 1, &y0[1], 1);
		RNP::TLASupport::LUSolve<'N'>(l-1, 1, &zz[1+1*(l+1)], l+1, iwork, &y0[1], l+1);
		y0[l] = 0.;

		yl[0] = 0.;
		RNP::TBLAS::Copy(l-1, &z[1+l*(l+1)], 1, &yl[1], 1);
		RNP::TLASupport::LUSolve<'N'>(l-1, 1, &zz[1+1*(l+1)], l+1, iwork, &yl[1], l+1);
		yl[l] = -1.;
		// Convex combination
		RNP::TBLAS::MultHermV<'U'>(l+1, 1., z, l+1, y0, 1, 0., y, 1);
		std::complex<double> kappa0 = sqrt(RNP::TBLAS::ConjugateDot(l+1, y0, 1, y, 1));

		RNP::TBLAS::MultHermV<'U'>(l+1, 1., z, l+1, yl, 1, 0., y, 1);
		std::complex<double> kappal = sqrt(RNP::TBLAS::ConjugateDot(l+1, yl, 1, y, 1));

		RNP::TBLAS::MultHermV<'U'>(l+1, 1., z, l+1, y0, 1, 0., y, 1);
		std::complex<double> varrho = RNP::TBLAS::ConjugateDot(l+1, yl, 1, y, 1) / (kappa0*kappal);

		double factor = std::abs(varrho); if(params.angle > factor){ factor = params.angle; }
		std::complex<double> hatgamma =  varrho/abs(varrho)*factor * (kappa0/kappal);

		RNP::TBLAS::Axpy(l+1, -hatgamma, yl, 1, y0, 1);
		// Update
		omega = y0[l];

		RNP::TBLAS::MultMV<'N'>(n, l, -1., &u[0+1*n], n, &y0[1], 1, 1., &u[0+0*n], 1);
		RNP::TBLAS::MultMV<'N'>(n, l,  1., &r[0+0*n], n, &y0[1], 1, 1., x, 1);
		RNP::TBLAS::MultMV<'N'>(n, l, -1., &r[0+1*n], n, &y0[1], 1, 1., r, 1);

		RNP::TBLAS::MultHermV<'U'>(l+1, 1., z, l+1, y0, 1, 0., y, 1);
		rnrm = std::sqrt(std::abs(RNP::TBLAS::ConjugateDot(l+1, y0, 1, y, 1)));
		// The reliable update part
		if(rnrm > mxnrmx){ mxnrmx = rnrm; }
		if(rnrm > mxnrmr){ mxnrmr = rnrm; }
		bool xpdt = (rnrm < params.delta*rnrm0 && rnrm0 < mxnrmx);
		if((rnrm < params.delta*mxnrmr && rnrm0 < mxnrmr) || xpdt){
			A(x, r, Adata);
			if(NULL != P){
				P(r, Pdata);
			}
			++nmv;
			for(size_t i = 0; i < n; ++i){
				r[i+0*n] =  bp[i+0*n] - r[i+0*n];
			}
			mxnrmr = rnrm;
			if(xpdt){
				RNP::TBLAS::Axpy(n, 1., x, 1, xp, 1);
				RNP::TBLAS::SetMatrix<'N'>(n, 1, 0., 0., x, 1);
				RNP::TBLAS::Copy(n, r, 1, bp, 1);
				mxnrmx = rnrm;
			}
		}
	}
	// End of iterations
	RNP::TBLAS::Axpy(n, 1., xp, 1, x, 1);
	// Check stopping criterion
	A(x, r, Adata);
	for(size_t i = 0; i < n; ++i){
		r[i+0*n] = b[i] - r[i+0*n];
	}
	if(NULL != P){
		P(r, Pdata);
	}
	rnrm = RNP::TBLAS::Norm2(n, r, 1);
	if (rnrm > ls_params.tol*rnrm0) info = 1;

	ls_params.max_iterations = iter;
	ls_params.tol = rnrm/rnrm0;
	ls_params.max_mult = nmv;

	if(NULL == workspace){
		RNP::BiCGSTABl_free(n, params, _work);
		delete _work;
	}
	return info;
}

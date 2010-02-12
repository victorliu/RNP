#include <cstddef>
#include <cmath>
#include <complex>
#include "TBLAS.h"
#include "TLASupport.h"
#include "LinearSolve.h"
#include "IterativeLinearSolvers.h"

// Implementation of GMRES - Generalized Minimum Residual

// According to the author, Gerard Sleijpen:
// The code is distributed under the terms of the GNU General Public
// License (version 2 of the License, or any later version) as
// published by the Free Software Foundation.
//
// (However, much of the code has been rewritten.)

void RNP::GMRES_alloc(size_t n, const GMRES_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	const size_t m = params.m;
	workspace->zwork = new std::complex<double>[(2+m+n)*(m+1)+m];
	workspace->rwork = new double[m];
}
void RNP::GMRES_free(size_t n, const GMRES_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	if(NULL != workspace->zwork){ delete [] workspace->zwork; workspace->zwork = NULL; }
	if(NULL != workspace->rwork){ delete [] workspace->rwork; workspace->rwork = NULL; }
}

int RNP::GMRES(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	RNP::LinearSolverParams &ls_params,
	const RNP::GMRES_params &params,
	// Optional parameters
	RNP::LinearSolverWorkspace *workspace,
	//std::complex<double> *work, // if not null, should be length (2+m+n)*(m+1)+m, m = params.m
	//double *rwork, // if not null, should be length params.m
	void *Adata,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*),
	void *Pdata
){
	const size_t mxm = params.m;
	
	// Set initial x
	if(!ls_params.x_initialized){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
	}
	if(NULL != P){ P(b, Pdata); }
	
	RNP::LinearSolverWorkspace *_work = workspace;
	if(NULL == workspace){
		_work = new RNP::LinearSolverWorkspace();
		RNP::GMRES_alloc(n, params, _work);
	}
	double *cs = _work->rwork;
	std::complex<double> *sn = _work->zwork; // length mxm
	std::complex<double> *rs = sn + mxm; // length mxm+1
	std::complex<double> *y = rs + mxm+1; // length mxm+1
	std::complex<double> *v = y + mxm+1; // length n*(mxm+1)
	std::complex<double> *hh = v + n*(mxm+1); // length mxm*(mxm+1)
	
	double rnrm0 = RNP::TBLAS::Norm2(n, b, 1);
	double rnrm = rnrm0;
	double eps1 = ls_params.tol * rnrm;
	RNP::TBLAS::Copy(n, b, 1, &v[0+0*n], 1);
	
	size_t iter = 0;
	size_t num_multmv = 0;
	
	size_t maxit = ls_params.max_iterations;
	if(0 == maxit){
		maxit = 2*n;
		if(1000 < maxit){ maxit = 1000; }
	}
	if(0 == ls_params.max_mult){ ls_params.max_mult = (size_t)-1; }
	
	// iteration
	while(num_multmv <= ls_params.max_mult && iter <= maxit && rnrm > eps1){
		RNP::TBLAS::Scale(n, 1./rnrm, &v[0+0*n], 1);
		rs[0] = rnrm;
		
		size_t m = 0;
		while(m < mxm && rnrm > eps1 && num_multmv <= ls_params.max_mult){
			A(&v[0+m*n], &v[0+(m+1)*n], Adata); ++num_multmv;
			if(NULL != P){ P(&v[0+(m+1)*n], Pdata); }
			for(size_t i = 0; i <= m; ++i){
				hh[i+m*(mxm+1)] = RNP::TBLAS::ConjugateDot(n, &v[0+i*n], 1, &v[0+(m+1)*n], 1);
				RNP::TBLAS::Axpy(n, -hh[i+m*(mxm+1)], &v[0+i*n], 1, &v[0+(m+1)*n], 1);
			}
			hh[m+1+m*(mxm+1)] = RNP::TBLAS::Norm2(n, &v[0+(m+1)*n], 1);
			RNP::TBLAS::Scale(n, 1./hh[m+1+m*(mxm+1)], &v[0+(m+1)*n], 1);
			for(size_t i = 0; i < m; ++i){
				RNP::TLASupport::ApplyPlaneRotation(1, &hh[i+m*(mxm+1)], 1, &hh[i+1+m*(mxm+1)], 1, cs[i], sn[i]);
			}
			std::complex<double> rcs;
			RNP::TLASupport::GeneratePlaneRotation(hh[m+m*(mxm+1)], hh[m+1+m*(mxm+1)], &cs[m], &sn[m], &rcs);
			hh[m+m*(mxm+1)] = rcs;
			hh[m+1+m*(mxm+1)] = 0;
			rs[m+1] = 0;
			RNP::TLASupport::ApplyPlaneRotation(1, &rs[m], 1, &rs[m+1], 1, cs[m], sn[m]);
			rnrm = std::abs(rs[m+1]);
			
			++m;
		}

		// compute approximate solution x
		RNP::TBLAS::Copy(m, rs, 1, y, 1);
		RNP::TBLAS::SolveTrV<'U','N','N'>(m, hh, mxm+1, y, 1);
		RNP::TBLAS::MultMV<'N'>(n, m, 1., v, n, y, 1, 1., x, 1);

		// compute residual for restart
		A(x, &v[0+1*n], Adata);
		if(NULL != P){ P(&v[0+1*n], Pdata); }
		RNP::TBLAS::Copy(n, b, 1, &v[0+0*n], 1);
		RNP::TBLAS::Axpy(n, -1., &v[0+1*n], 1, &v[0+0*n], 1);
		rnrm = RNP::TBLAS::Norm2(n, &v[0+0*n], 1);
		++iter;
	}

	ls_params.max_iterations = iter;
	ls_params.tol = rnrm / rnrm0;
	ls_params.max_mult = num_multmv;
	
	if(NULL == workspace){
		RNP::GMRES_free(n, params, _work);
		delete _work;
	}
	return 0;
}

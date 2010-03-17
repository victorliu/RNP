#include <cstddef>
#include <cmath>
#include <complex>
#include <cstdlib>
#include "TBLAS.h"
#include "TLASupport.h"
#include "LinearSolve.h"
#include "IterativeLinearSolvers.h"

// Implementation of IDR(s) - Induced Dimension Reduction method

// Original copyright:
// The software is distributed without any warranty.
// Martin van Gijzen and Peter Sonneveld
// Copyright (c) December 2008
//
// (Note that many changes have been made to the original code.)

void RNP::IDRs_alloc(size_t n, const IDRs_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	const size_t s = params.s;
	workspace->zwork = new std::complex<double>[(3*n + 2*s)*(s+1)];
	workspace->iwork = new size_t[s];
}
void RNP::IDRs_free(size_t n, const IDRs_params &params, LinearSolverWorkspace *workspace){
	if(NULL == workspace){ return; }
	if(NULL != workspace->zwork){ delete [] workspace->zwork; workspace->zwork = NULL; }
	if(NULL != workspace->iwork){ delete [] workspace->iwork; workspace->iwork = NULL; }
}

int RNP::IDRs(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	RNP::LinearSolverParams &ls_params,
	const RNP::IDRs_params &params,
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
	const double normb = RNP::TBLAS::Norm2(n, b, 1);
	if(0 == normb){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
		return 0;
	}
	const double tolb = ls_params.tol*normb; // compute tolerance
	size_t num_multmv = 0;
	
	// Set initial x
	if(!ls_params.x_initialized){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
	}
	
	const size_t s = params.s;
	
	RNP::LinearSolverWorkspace *_work = workspace;
	if(NULL == workspace){
		_work = new RNP::LinearSolverWorkspace();
		RNP::IDRs_alloc(n, params, _work);
	}
	size_t *ipiv = _work->iwork;
	
	std::complex<double> *r = _work->zwork; // length n
	A(x,r,Adata); ++num_multmv;
	for(size_t i = 0; i < n; ++i){ r[i] = b[i]-r[i]; }
	double normr = RNP::TBLAS::Norm2(n, r, 1);
	// Now, r = b-A*x
	
	std::complex<double> *Q = r+n; // length n*s;
	{ // set up shadow space
		
		for(size_t j = 0; j < s; ++j){
			for(size_t i = 0; i < n; ++i){
				Q[i+j*n] = (double)rand()/(double)RAND_MAX - 0.5;
			}
		}
		// Orthogonalize Q
		for(size_t j = 0; j < s; ++j){
			double inorm = double(1)/RNP::TBLAS::Norm2(n, &Q[0+j*n], 1);
			for(size_t i = 0; i < n; ++i){
				Q[i+j*n] *= inorm;
			}
			for(size_t k = j+1; k < s; ++k){
				std::complex<double> dot = 0;
				for(size_t i = 0; i < n; ++i){
					dot += std::conj(Q[i+j*n]) * Q[i+k*n];
				}
				for(size_t i = 0; i < n; ++i){
					Q[i+k*n] -= dot*Q[i+j*n];
				}
			}
			//std::cout << "Q[" << j << "] = "; print(n, &Q[0+j*n]); std::cout << std::endl;
		}
	}
	
	std::complex<double> *G = Q+n*s; // length n*s
	std::complex<double> *U = G+n*s; // length n*s
	std::complex<double> *M = U+n*s; // length s*s
	std::complex<double> *Mcopy = M+s*s; // length s*s

	for(size_t j = 0; j < s; ++j){
		for(size_t i = 0; i < n; ++i){
			G[i+j*n] = 0;
			U[i+j*n] = 0;
		}
		for(size_t i = 0; i < s; ++i){
			if(i == j){
				M[i+j*s] = 1;
			}else{
				M[i+j*s] = 0;
			}
		}
	}
	std::complex<double> *f = Mcopy+s*s; // length s
	std::complex<double> *c = f+s; // length s
	std::complex<double> *v = c+s; // length n
	std::complex<double> *t = v+n; // length n
	size_t iter = 0;
	std::complex<double> om = 1;
	
	size_t maxit = ls_params.max_iterations;
	if(0 == maxit){
		maxit = 2*n;
		if(1000 < maxit){ maxit = 1000; }
	}
	if(0 == ls_params.max_mult){ ls_params.max_mult = (size_t)-1; }
	
	int ret = 0;
	while(normr > tolb && iter < maxit && num_multmv <= ls_params.max_mult){
//		std::cout << "iter = " << iter << std::endl;
		
		// generate RHS for small system
		for(size_t j = 0; j < s; ++j){
			std::complex<double> sum = 0;
			for(size_t i = 0; i < n; ++i){
				sum += r[i] * std::conj(Q[i+j*n]);
			}
			f[j] = sum;
		}
		
		for(size_t k = 0; k < s; ++k){
			// solve small systems of M(k:s,k:s)*c(k:s) = f(k:s)
			{
				// Copy over stuff for a destructive LU solve in Mcopy
				for(size_t j = k; j < s; ++j){
					for(size_t i = k; i < s; ++i){
						Mcopy[i+j*s] = M[i+j*s];
					}
					c[j] = f[j];
				}
				// Perform LU solve...
				RNP::LinearSolve<'N'>(s-k, 1, &Mcopy[k+k*s], s, &c[k], s, NULL, ipiv);
			}
			// v = r - G(:,k:s)*c;
			for(size_t i = 0; i < n; ++i){
				std::complex<double> sum = 0;
				for(size_t j = k; j < s; ++j){
					sum += G[i+j*n]*c[j];
				}
				v[i] = r[i] - sum;
			}
			if(NULL != P){
				P(v, Pdata);
			}
			
			//U(:,k) = U(:,k:s)*c + om*v;
			for(size_t i = 0; i < n; ++i){
				std::complex<double> sum = 0;
				for(size_t j = k; j < s; ++j){
					sum += U[i+j*n]*c[j];
				}
				U[i+k*n] = sum + om*v[i];
			}
			//G(:,k) = A*U(:,k);
			A(&U[0+k*n], &G[0+k*n], Adata); ++num_multmv;
			
			// Bi-Orthogonalise the new basis vectors
			for(size_t j = 0; j < k; ++j){
				std::complex<double> alpha = 0;
				for(size_t i = 0; i < n; ++i){
					alpha += std::conj(Q[i+j*n])*G[i+k*n];
				}
				alpha /= M[j+j*s];
				for(size_t i = 0; i < n; ++i){
					G[i+k*n] -= alpha*G[i+j*n];
				}
				for(size_t i = 0; i < n; ++i){
					U[i+k*n] -= alpha*U[i+j*n];
				}
			}
			// New column of M = (Q'*G)'  (first k-1 entries are zero)
			for(size_t j = k; j < s; ++j){
				std::complex<double> sum = 0;
				for(size_t i = 0; i < n; ++i){
					sum += G[i+k*n]*std::conj(Q[i+j*n]);
				}
				M[j+k*s] = sum;
			}

			// Make r orthogonal to p_i, i = 1..k
			std::complex<double> beta = f[k]/M[k+k*s];
			for(size_t i = 0; i < n; ++i){
				r[i] -= beta*G[i+k*n];
			}
			for(size_t i = 0; i < n; ++i){
				x[i] += beta*U[i+k*n];
			}

			++iter;
			normr = RNP::TBLAS::Norm2(n, r, 1);

			if(normr < tolb || iter == maxit){ break; }
			
			// New f = Q'*r (first k  components are zero)
			for(size_t j = k+1; j < s; ++j){
				f[j] -= beta*M[j+k*s];
			}
		} // end k loop
		
		// If we break'd out of the inner loop, do so again
		if(normr < tolb){ break; }

		// Now we have sufficient vectors in G_j to compute residual in G_j+1
		// Note: r is already perpendicular to Q so v = r
		for(size_t i = 0; i < n; ++i){ v[i] = r[i]; }
		if(NULL != P){
			P(v, Pdata);
		}
		A(v, t, Adata); ++num_multmv;
		{ // compute new omega
			double norms = RNP::TBLAS::Norm2(n, r, 1), normt = RNP::TBLAS::Norm2(n, t, 1);
			std::complex<double> ts = 0;
			for(size_t i = 0; i < n; ++i){
				ts += std::conj(t[i])*r[i];
			}
			double rho = std::abs(ts/(normt*norms));
			om = ts/(normt*normt);
			if(rho < params.angle){
				om *= params.angle/rho;
			}
		}
		
		for(size_t i = 0; i < n; ++i){ r[i] -= om*t[i]; }
		for(size_t i = 0; i < n; ++i){ x[i] += om*v[i]; }
		normr = RNP::TBLAS::Norm2(n, r, 1);
		++iter;
	}
	
	ls_params.max_iterations = iter;
	ls_params.tol = normr/normb;
	ls_params.max_mult = num_multmv;
	
	if(NULL == workspace){
		RNP::IDRs_free(n, params, _work);
		delete _work;
	}
	return ret;
}

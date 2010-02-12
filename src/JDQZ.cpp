#include <cmath>
#include <complex>
#include <limits>
#include <ctime>
#include "TBLAS.h"
#include "Random.h"
#include "TLASupport.h"
#include "GeneralizedEigensystems.h"
#include "IterativeGeneralizedEigensystems.h"
#include "IterativeLinearSolvers.h"

// According to the author of JDQZ, Gerard Sleijpen:
// The code is distributed under the terms of the GNU General Public
// License (version 2 of the License, or any later version) as
// published by the Free Software Foundation.
//
// (However, much of the code has been rewritten.)

bool RNP::JDQZSort_MaxReal(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	// If n1 = a+bi, d1 = c+di,
	// n1/d1 = [(ac+bd) + i(bc-ad)] / (cc+dd)

	// naive method
	return (n1/d1).real() > (n2/d2).real();
}
bool RNP::JDQZSort_MaxAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	return std::abs(n1/d1) > std::abs(n2/d2);
}
bool RNP::JDQZSort_MinReal(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	return (n1/d1).real() > (n2/d2).real();
}
bool RNP::JDQZSort_MinAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	return std::abs(n1/d1) < std::abs(n2/d2);
}
bool RNP::JDQZSort_MaxImag(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	return (n1/d1).imag() > (n2/d2).imag();
}
bool RNP::JDQZSort_MinImag(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	return (n1/d1).imag() < (n2/d2).imag();
}
bool RNP::JDQZSort_Nearest(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	const std::complex<double> *target = reinterpret_cast<std::complex<double>*>(data);
	std::complex<double> z1 = n1/d1 - *target;
	std::complex<double> z2 = n2/d2 - *target;
	return std::abs(z1) < std::abs(z2);
}
bool RNP::JDQZSort_NearestAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data){
	const double *target = reinterpret_cast<double*>(data);
	double a1 = std::abs(n1/d1) - *target;
	double a2 = std::abs(n2/d2) - *target;
	return std::abs(a1) < std::abs(a2);
}


// Sorts a QZ decomposition
// Inputs:
//   A generalized Schur decomposition of (A,B) such that
//   U*S*V^H == A and U*T*V^H == B
//   where k is the dimension of the problem.
//
//   EvalSorter takes (num1,den1, num2,den2) and returns
//   true if num1/den1 should come before num2/den2.
// Return:
//   The eigenvalues are sorted according to the comparison
//   function EvalSorter. This must be a strict weak ordering.
void QZSort(size_t k,
	std::complex<double> *s, size_t lds,
	std::complex<double> *t, size_t ldt,
	std::complex<double> *u, size_t ldu,
	std::complex<double> *v, size_t ldv,
	bool (*EvalSorter)(const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, void*),
	void *sorter_data
){
	for(size_t i = 0; i < k; ++i){
		size_t j = i;
		// Determine the order of the diagonal elements from i:k-1 (0-based inclusive)
		for(size_t l = i+1; l < k; ++l){
			if(EvalSorter(s[l+l*lds], t[l+l*ldt], s[j+j*lds], t[j+j*ldt], sorter_data)){
				j = l;
			}
		}
//		for(size_t qq = 0; qq < k; ++qq){
//			std::cout << (s[qq+qq*lds]/t[qq+qq*ldt]).real() << std::endl;
//		}std::cout << std::endl;
		if(j != i){
//			std::cout << "Swapping " << i << " with " << j << std::endl;
			// Perform a swap between i and j; we always have i < j
			//for(size_t p = i; p < j; ++p){
			for(int p = j-1; p >= (int)i; --p){
				double f = std::max(std::abs(t[p+1+(p+1)*ldt]), std::abs(s[p+1+(p+1)*lds]));
				std::complex<double> ft = t[p+1+(p+1)*ldt] / f;
				std::complex<double> fs = s[p+1+(p+1)*lds] / f;
				bool tlts = (std::abs(ft) < std::abs(fs));
				std::complex<double> c1 = ft*s[p+p*lds] - fs*t[p+p*ldt];
				std::complex<double> c2 = ft*s[p+(p+1)*lds] - fs*t[p+(p+1)*ldt];

				double cs;
				std::complex<double> sn, r;
				RNP::TLASupport::GeneratePlaneRotation(c2, c1, &cs, &sn, &r);
				RNP::TLASupport::ApplyPlaneRotation(p+2, &s[0+(p+1)*lds], 1, &s[0+p*lds], 1, cs, sn);
				RNP::TLASupport::ApplyPlaneRotation(p+2, &t[0+(p+1)*ldt], 1, &t[0+p*ldt], 1, cs, sn);
				RNP::TLASupport::ApplyPlaneRotation(k, &v[0+(p+1)*ldu], 1, &v[0+p*ldu], 1, cs, sn);
				if(tlts){
					c1 = s[p+p*lds];
					c2 = s[p+1+p*lds];
				}else{
					c1 = t[p+p*ldt];
					c2 = t[p+1+p*ldt];
				}
				RNP::TLASupport::GeneratePlaneRotation(c1, c2, &cs, &sn, &r);
				RNP::TLASupport::ApplyPlaneRotation(k-p, &s[p+p*lds], lds, &s[p+1+p*lds], lds, cs, sn);
				RNP::TLASupport::ApplyPlaneRotation(k-p, &t[p+p*ldt], ldt, &t[p+1+p*ldt], ldt, cs, sn);
				RNP::TLASupport::ApplyPlaneRotation(k, &u[0+p*ldv], 1, &u[0+(p+1)*ldv], 1, cs, std::conj(sn));
			}
		}
	}
}

// Input:
//   v is an nxk matrix (ldv = n) with orthonormal columns
//   w is a vector of length n
// Output:
//   w is orthogonalized against each column of v, and optionally normalized
int ModifiedGramSchmidt(size_t n, size_t k, const std::complex<double> *v, std::complex<double> *w, bool normalize_w = true){
	double s1 = RNP::TBLAS::Norm2(n, w, 1);
	int info = 0;
	for(size_t i = 0; i < k; ++i){
		double s0 = s1;
		std::complex<double> znrm = RNP::TBLAS::ConjugateDot(n, &v[0+i*n], 1, w, 1);
		RNP::TBLAS::Axpy(n, -znrm, &v[0+i*n], 1, w, 1);
		s1 = RNP::TBLAS::Norm2(n, w, 1);
		if(s1 <= s0/100){
			// If the norm of w was reduced drastically,
			// w was almost parallel to a vector in v
			s0 = s1;
			std::complex<double> ztmp = RNP::TBLAS::ConjugateDot(n, &v[0+i*n], 1, w, 1);
			znrm += ztmp;
			RNP::TBLAS::Axpy(n, -ztmp, &v[0+i*n], 1, w, 1);
			s1 = RNP::TBLAS::Norm2(n, w, 1);
			if(s1 <= s0/100){
				// we have a problem...
				if(0 == info){
					info = 1 + i;
				}
			}
		}
	}
	if(normalize_w){
		RNP::TBLAS::Scale(n, 1./s1, w, 1);
	}
	return info;
}

struct jdqz_op_data{
	size_t n;
	void (*Aop)(const std::complex<double> *, std::complex<double> *, void*);
	void (*Bop)(const std::complex<double> *, std::complex<double> *, void*);
	std::complex<double> *work;
	std::complex<double> alpha, beta;
	void *Adata, *Bdata;
};
// Computes y = beta*Ax - alpha*Bx
void jdqz_op(const std::complex<double> *x, std::complex<double> *y, void *_data){
	jdqz_op_data *data = reinterpret_cast<jdqz_op_data*>(_data);
	if(&RNP::JDQZ_Identity == data->Bop){
		// B is identity, so we only need to comput y = beta*Ax
		data->Aop(x, y, data->Adata);
		RNP::TBLAS::Scale(data->n, data->beta, y, 1);
	}else{
		data->Aop(x, data->work, data->Adata);
		data->Bop(x, y, data->Bdata);
		for(size_t i = 0; i < data->n; ++i){
			y[i] = data->beta * data->work[i] - data->alpha * y[i];
		}
	}
}

void RNP::JDQZ_Identity(const std::complex<double> *x, std::complex<double> *y, void *_data){
	size_t *n = reinterpret_cast<size_t*>(_data);
	RNP::TBLAS::Copy(*n, x, 1, y, 1);
}

struct jdqz_preconditioner_data{
	size_t n;
	void (*preconditioner)(std::complex<double> *, void *);
	void *preconditioner_data;
	size_t nq; // number of columns in q
	std::complex<double> *q;
	std::complex<double> *kz;
	std::complex<double> *invqkz;
	size_t ldqkz; // typically jmax
	size_t *ipiv; // length ldqkz; should store pivots for LU of invqkz
	std::complex<double> *work; // length n
};

void jdqz_preconditioner(std::complex<double> *x, void *_data){
	jdqz_preconditioner_data *data = reinterpret_cast<jdqz_preconditioner_data*>(_data);
	if(NULL != data->preconditioner){
		data->preconditioner(x, data->preconditioner_data);
	}
	RNP::TBLAS::MultMV<'C'>(data->n, data->nq, 1., data->q, data->n, x, 1, 0., data->work, 1);
	RNP::TLASupport::LUSolve<'N'>(data->nq, 1, data->invqkz, data->ldqkz, data->ipiv, data->work, data->ldqkz);
	RNP::TBLAS::MultMV<'N'>(data->n, data->nq, -1., data->kz, data->n, data->work, 1, 1., x, 1);
}

struct jdqz_linsolve_data{
	RNP::LinearSolverWorkspace *workspace;
	RNP::LinearSolverParams ls_params;
	jdqz_linsolve_data(){
		workspace = NULL;
	}
};

void RNP::JDQZ_IDRs_adapter(
	size_t n,
	void (*Op)(const std::complex<double> *, std::complex<double> *, void*),
	void *op_data,
	std::complex<double> *x,
	std::complex<double> *b,
	double tol,
	size_t &max_mult,
	void (*jdqz_precon)(std::complex<double> *x, void *precon_data),
	void *jdqz_precon_data,
	void *_data
){
	jdqz_linsolve_data *data = reinterpret_cast<jdqz_linsolve_data*>(_data);
	static RNP::IDRs_params params;
	data->ls_params.max_mult = max_mult;
	data->ls_params.tol = tol;
	RNP::TBLAS::Scale(n, -1., b, 1);
	RNP::IDRs(n, Op, b, x, data->ls_params, params, data->workspace, op_data, jdqz_precon, jdqz_precon_data);
	max_mult = data->ls_params.max_mult;
}

void RNP::JDQZ_GMRES_adapter(
	size_t n,
	void (*Op)(const std::complex<double> *, std::complex<double> *, void*),
	void *op_data,
	std::complex<double> *x,
	std::complex<double> *b,
	double tol,
	size_t &max_mult,
	void (*jdqz_precon)(std::complex<double> *x, void *precon_data),
	void *jdqz_precon_data,
	void *_data
){
	jdqz_linsolve_data *data = reinterpret_cast<jdqz_linsolve_data*>(_data);
	static RNP::GMRES_params params;
	data->ls_params.max_mult = max_mult;
	data->ls_params.tol = tol;
	RNP::TBLAS::Scale(n, -1., b, 1);
	RNP::GMRES(n, Op, b, x, data->ls_params, params, data->workspace, op_data, jdqz_precon, jdqz_precon_data);
	max_mult = data->ls_params.max_mult;
}

void jdqz_pick_search_space(size_t k, size_t &jmin, size_t &jmax){
	if(0 == jmin){ jmin = 2*k; }
	if(0 == jmax){
		jmax = 3*k;
		if(jmax < 20){
			jmax = 20;
		}
	}
}

// Inputs:
//   n - size of the problem
//   k - number of eigenstuffs to compute
//   Aop - A routine that applies the A matrix to the first argument and stores it in the second
//   Bop - just like Aop
//   EvalSorter - sorter for eigenvalues (args: alpha1,beta1,alpha2,beta2,data)
//   preconditioner - Applies (approximatino of) inv(A-tau*B) to its argument
// Ouptuts:
//   alpha,beta - numerator/denominator pairs for eigenvalues (length k)
//   vr,ldvr - eigenvalue matrix and leading dimension
//
int RNP::JDQZ(
	size_t n, size_t kmax,
	void (*Aop)(const std::complex<double> *, std::complex<double> *, void*),
	void (*Bop)(const std::complex<double> *, std::complex<double> *, void*),
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *eivec, size_t ldeivec,
	const std::complex<double> &target,
	bool (*EvalSorter)(const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, const std::complex<double> &, void*),
	void (*preconditioner)(std::complex<double> *, void *),
	const RNP::JDQZParams &params,
	RNP::JDQZWorkspace *workspace,
	// User data parameters passed onto the corresponding arguments
	void *Aop_data,
	void *Bop_data,
	void *EvalSorter_data,
	void *preconditioner_data
){
	// Check inputs
	if(kmax > n){ return -2; }
	if(ldeivec < n){ return -8; }
	size_t jmin = params.min_search_space;
	size_t jmax = params.max_search_space;
	jdqz_pick_search_space(kmax, jmin, jmax);
	if(jmax <= jmin){ return -11; }
	if(kmax > jmax){ return -11; }

	// Setup the workspace
	RNP::JDQZWorkspace *_work = workspace;
	if(NULL == workspace){
		_work = new RNP::JDQZWorkspace();
		RNP::JDQZ_alloc(n, kmax, params, _work);
	}
	std::complex<double> *rhs, *temp, *v, *w, *av, *bv, *aux, *q, *z, *kz;
	std::complex<double> *f, *ma, *mb, *ra, *rb, *zma, *zmb, *vsl, *vsr;
	std::complex<double> *mqkz, *invqkz, *aconv, *bconv;
	std::complex<double> *zwork; // for zggev

	// leading dimension n:
	//   rhs = zeros(n,1) // RHS of solves
	//   temp = zeros(n,1) // workspace for computing (alpha*A-beta*B)*x, general temporary storage for a vector
	//   //DO-NOT-COUNT u = linear_solver_space (for GMRES, zeros(n,m))
	//   v = zeros(n,jmax) // pointer to search space JDQZ
	//   w = zeros(n,jmax) // pointer to test subspace JDQZ
	//   av = zeros(n,jmax) // pointer to subspace AV
	//   bv = zeros(n,jmax) // pointer to subspace BV
	//   aux = zeros(n,jmax) // rotating buffer space; points at v,w,av,bv in rotation to prevent unnecessary copies
	//   q = zeros(n,kmax) // pointer to search Schur basis in JDQZ
	//   z = zeros(n,kmax) // pointer to test Schur basis in JDQZ
	//   kz = zeros(n,kmax) // stores matrix inv(K)*Z_k
	// leading dimension jmax:
	//   f = zeros(jmax,1) // temporary space used in the preconditioner
	//   ma = zeros(jmax,jmax) // M_A matrix (pre-Schur decomposition matrix)
	//   mb = zeros(jmax,jmax) // M_B matrix
	//   ra = zeros(jmax,jmax) // r_A matrix (stores up residuals)
	//   rb = zeros(jmax,jmax) // r_B matrix
	//   zma = zeros(jmax,jmax) // M_A matrix (post-Schur decomposition matrix)
	//   zmb = zeros(jmax,jmax) // M_B matrix
	//   vsl = zeros(jmax,jmax) // left Schur vectors
	//   vsr = zeros(jmax,jmax) // right Schur vectors
	//   mqkz = zeros(jmax,jmax)
	//   invqkz = zeros(jmax,jmax)
	//   aconv = zeros(jmax,1) // temp holding of eigenvalues (may overflow alpha/beta)
	//   bconv = zeros(jmax,1)
	//_work = new std::complex<double>[n*(2+5*jmax+3*kmax)+jmax*(10*jmax+5)];
	
	rhs = _work->zwork;
	temp = rhs + n;
	v = temp + n;
	w = v + n*jmax;
	av = w + n*jmax;
	bv = av + n*jmax;
	aux = bv + n*jmax;
	q = aux + n*jmax;
	z = q + n*kmax;
	kz = z + n*kmax;
	f = kz + n*kmax;
	ma = f + jmax;
	mb = ma + jmax*jmax;
	ra = mb + jmax*jmax;
	rb = ra + jmax*jmax;
	zma = rb + jmax*jmax;
	zmb = zma + jmax*jmax;
	vsl = zmb + jmax*jmax;
	vsr = vsl + jmax*jmax;
	mqkz = vsr + jmax*jmax;
	invqkz = mqkz + jmax*jmax;
	aconv = invqkz + jmax*jmax;
	bconv = aconv + jmax;
	zwork = bconv + jmax;
	
	/*
	rhs = new std::complex<double>[n];
	temp = new std::complex<double>[n];
	v = new std::complex<double>[n*jmax];
	w = new std::complex<double>[n*jmax];
	av = new std::complex<double>[n*jmax];
	bv = new std::complex<double>[n*jmax];
	aux = new std::complex<double>[n*jmax];
	q = new std::complex<double>[n*kmax];
	z = new std::complex<double>[n*kmax];
	kz = new std::complex<double>[n*kmax];
	f = new std::complex<double>[jmax];
	ma = new std::complex<double>[jmax*jmax];
	mb = new std::complex<double>[jmax*jmax];
	ra = new std::complex<double>[jmax*jmax];
	rb = new std::complex<double>[jmax*jmax];
	zma = new std::complex<double>[jmax*jmax];
	zmb = new std::complex<double>[jmax*jmax];
	vsl = new std::complex<double>[jmax*jmax];
	vsr = new std::complex<double>[jmax*jmax];
	mqkz = new std::complex<double>[jmax*jmax];
	invqkz = new std::complex<double>[jmax*jmax];
	aconv = new std::complex<double>[jmax];
	bconv = new std::complex<double>[jmax];
	*/
	if(NULL == Bop){
		Bop = &RNP::JDQZ_Identity;
		Bop_data = (void*)&n;
	}
	jdqz_op_data op_data;
	op_data.work = temp;
	op_data.n = n;
	op_data.Aop = Aop;
	op_data.Bop = Bop;
	op_data.Adata = Aop_data;
	op_data.Bdata = Bop_data;
	// rwork needed for zgges
	double *rwork = _work->rwork;
	size_t *ipivqkz = _work->iwork;

	jdqz_preconditioner_data precon_data;
	precon_data.n = n;
	precon_data.preconditioner = preconditioner;
	precon_data.preconditioner_data = preconditioner_data;
	precon_data.q = q;
	precon_data.kz = kz;
	precon_data.invqkz = invqkz;
	precon_data.ldqkz = jmax;
	precon_data.ipiv = ipivqkz;
	precon_data.work = f;

	size_t num_multmv = 0;
	bool ok = true;
	double evcond = sqrt(std::norm(target) + 1);
	std::complex<double> shifta = target / evcond;
	std::complex<double> shiftb = 1. / evcond;
	std::complex<double> targeta = shifta;
	std::complex<double> targetb = shiftb;
	std::complex<double> zalpha = shifta;
	std::complex<double> zbeta = shiftb;
	size_t step = 0;
	int solvestep = 0;
	size_t j = 0;
	size_t k = 0;

	int iseed[4] = { 3,3,1966,29 };

	while(k < kmax && step < params.max_iters){
		++step;
		++solvestep;
		if(j == 0){ // set v to initial value v0
			for(size_t i = 0; i < n; ++i){
				double ra[2];
				RNP::Random::RandomRealsUniform01(2, ra, iseed);
				v[i+0*n] = 2*ra[0]-1;
			}
			for(size_t i = 0; i < n; ++i){
				double ra[2];
				RNP::Random::RandomRealsUniform01(2, ra, iseed);
				w[i+0*n] = 2*ra[0]-1;
			}
		}else{
			size_t linsolve_max_mult = params.max_mult;
			double deps = pow(2.0, -solvestep); // stopping tolerance for the solve
			precon_data.nq = k+1;
			if(j < jmin){
				linsolve_max_mult = 1;
			}
			op_data.alpha = zalpha;
			op_data.beta = zbeta;
			params.linear_solver(n, &jdqz_op, (void*)&op_data, &v[0+j*n], rhs, deps, linsolve_max_mult, &jdqz_preconditioner, (void*)&precon_data, params.linear_solver_data);
		}
		
		// The projected problem
		ModifiedGramSchmidt(n, j, v, &v[0+j*n]); // v = mgs(V,v); v = v/|v|
		ModifiedGramSchmidt(n, k, q, &v[0+j*n]); // this step was not mentioned in the paper
		if(params.testspace == 1){
			// Standard Petrov
			// testspace = alpha'*Av + beta'*Bv
			op_data.alpha = -std::conj(shiftb);
			op_data.beta = std::conj(shifta);
		}else if(params.testspace == 2){
			// Standard variable Petrov
			op_data.alpha = -std::conj(zbeta);
			op_data.beta = std::conj(zalpha);
		}else if(params.testspace == 3){
			// Harmonic Petrov
			// testspace = beta*Av - alpha*Bv
			op_data.alpha = shifta;
			op_data.beta = shiftb;
		}else if(params.testspace == 5){
			// testspace = Bv
			op_data.alpha = -1;
			op_data.beta = 0;
		}else if(params.testspace == 6){
			// testspace = Av
			op_data.alpha = 0;
			op_data.beta = 1;
		}
		jdqz_op(&v[0+j*n], &w[0+j*n], &op_data);
		
		ModifiedGramSchmidt(n, j, w, &w[0+j*n]); // w = mgs(W,w); w = w/|w|
		ModifiedGramSchmidt(n, k, z, &w[0+j*n]); // w = mgs(Z,w);
		Aop(&v[0+j*n], &av[0+j*n], Aop_data); // vA = A*v
		Bop(&v[0+j*n], &bv[0+j*n], Bop_data); // vB = B*v
		++j;
		// Make MA = [MA,W^* * vA; w^* * VA, w^* * vA]
		for(size_t jj = 0; jj < j; ++jj){
			for(size_t ii = 0; ii < j; ++ii){
				if(ii == j-1 || jj == j-1){
					ma[ii+jj*jmax] = RNP::TBLAS::ConjugateDot(n, &w[0+ii*n], 1, &av[0+jj*n], 1);
				}
				zma[ii+jj*jmax] = ma[ii+jj*jmax];
			}
		}
		for(size_t jj = 0; jj < j; ++jj){
			for(size_t ii = 0; ii < j; ++ii){
				if(ii == j-1 || jj == j-1){
					mb[ii+jj*jmax] = RNP::TBLAS::ConjugateDot(n, &w[0+ii*n], 1, &bv[0+jj*n], 1);
				}
				zmb[ii+jj*jmax] = mb[ii+jj*jmax];
			}
		}
		//RNP::GeneralizedSchurDecomposition(j, zma, jmax, zmb, jmax, alpha, beta, vsl, jmax, vsr, jmax, zwork, rwork);
		RNP::GeneralizedSchurDecomposition(j, zma, jmax, zmb, jmax, aconv, bconv, vsl, jmax, vsr, jmax, zwork, rwork);
		
		bool found = true;
		while(found){
			// Sort the Petrov pairs
			QZSort(j, zma, jmax, zmb, jmax, vsl, jmax, vsr, jmax, EvalSorter, EvalSorter_data);
			zalpha = zma[0];
			zbeta = zmb[0];
			evcond = sqrt(std::norm(zalpha) + std::norm(zbeta));

			// compute new q
			RNP::TBLAS::MultMV<'N'>(n, j, 1., v, n, vsr, 1, 0., &q[0+k*n], 1);
			ModifiedGramSchmidt(n, k, q, &q[0+k*n]);

			// compute new z
			RNP::TBLAS::MultMV<'N'>(n, j, 1., w, n, vsl, 1, 0., &z[0+k*n], 1);
			ModifiedGramSchmidt(n, k, z, &z[0+k*n]);

			// Make new qkz
			RNP::TBLAS::Copy(n, &z[0+k*n], 1, &kz[0+k*n], 1);
			if(NULL != preconditioner){
				preconditioner(&kz[0+k*n], preconditioner_data);
			}
			for(size_t jj = 0; jj <= k; ++jj){
				for(size_t ii = 0; ii <= k; ++ii){
					if(ii == k || jj == k){
						mqkz[ii+jj*jmax] = RNP::TBLAS::ConjugateDot(n, &q[0+ii*n], 1, &kz[0+jj*n], 1);
					}
					invqkz[ii+jj*jmax] = mqkz[ii+jj*jmax];
				}
			}
			RNP::TLASupport::LUDecomposition(k+1, k+1, invqkz, jmax, ipivqkz);

			// compute new (right) residual= beta Aq - alpha Bq and orthogonalize this vector on Z.
			op_data.alpha = zalpha;
			op_data.beta = zbeta;
			jdqz_op(&q[0+k*n], rhs, &op_data);
			ModifiedGramSchmidt(n, k, z, rhs, false);
			double rnrm = RNP::TBLAS::Norm2(n, rhs, 1) / evcond;
			if(rnrm < params.lock && ok){
				targeta = zalpha;
				targetb = zbeta;
				ok = false;
			}
			if(found = (rnrm < params.eps && (j > 1 || k == kmax - 1))){
				// store the eigenvalue
				alpha[k] = zalpha;
				beta[k] = zbeta;
				// increase the number of found evs by 1
				++k;
				solvestep = 0;
				if(k == kmax){
					break;
				}
				RNP::TBLAS::MultMM<'N','N'>(n, j-1, j, 1., v , n, &vsr[0+1*jmax], jmax, 0., aux, n);
				std::swap(v, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, j-1, j, 1., av, n, &vsr[0+1*jmax], jmax, 0., aux, n);
				std::swap(av, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, j-1, j, 1., bv, n, &vsr[0+1*jmax], jmax, 0., aux, n);
				std::swap(bv, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, j-1, j, 1., w , n, &vsl[0+1*jmax], jmax, 0., aux, n);
				std::swap(w, aux);
				--j;
				RNP::TBLAS::CopyMatrix<'A'>(j, j, &zma[1+1*jmax], jmax, ma, jmax);
				RNP::TBLAS::CopyMatrix<'A'>(j, j, ma, jmax, zma, jmax);
				RNP::TBLAS::CopyMatrix<'A'>(j, j, &zmb[1+1*jmax], jmax, mb, jmax);
				RNP::TBLAS::CopyMatrix<'A'>(j, j, mb, jmax, zmb, jmax);
				RNP::TBLAS::SetMatrix<'A'>(j, j, 0., 1., vsr, jmax);
				RNP::TBLAS::SetMatrix<'A'>(j, j, 0., 1., vsl, jmax);
				targeta = shifta;
				targetb = shiftb;
				ok = true;
			}else if(j == jmax){
				RNP::TBLAS::MultMM<'N','N'>(n, jmin, j, 1., v , n, vsr, jmax, 0., aux, n);
				std::swap(v, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, jmin, j, 1., av, n, vsr, jmax, 0., aux, n);
				std::swap(av, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, jmin, j, 1., bv, n, vsr, jmax, 0., aux, n);
				std::swap(bv, aux);
				RNP::TBLAS::MultMM<'N','N'>(n, jmin, j, 1., w , n, vsl, jmax, 0., aux, n);
				std::swap(w, aux);
				j = jmin;
				RNP::TBLAS::CopyMatrix<'A'>(j, j, zma, jmax, ma, jmax);
				RNP::TBLAS::CopyMatrix<'A'>(j, j, zmb, jmax, mb, jmax);
				RNP::TBLAS::SetMatrix<'A'>(j, j, 0., 1., vsr, jmax);
				RNP::TBLAS::SetMatrix<'A'>(j, j, 0., 1., vsl, jmax);
			}
		}
	}

	// Did enough eigenpairs converge?
	kmax = k;
	if(NULL != eivec){
		// Compute the Schur matrices if the eigenvectors are
		// wanted, work(1,temp) is used for temporary storage
		// Compute RA:
		RNP::TBLAS::SetMatrix<'L'>(k, k, 0., 0., ra, jmax);
		for(size_t i = 0; i < k; ++i){
			Aop(&q[0+i*n], temp, Aop_data);
			RNP::TBLAS::MultMV<'C'>(n, i+1, 1., z, n, temp, 1, 0., &ra[0+i*jmax], 1);
		}
		// Compute RB:
		RNP::TBLAS::SetMatrix<'L'>(k, k, 0., 0., rb, jmax);
		for(size_t i = 0; i < k; ++i){
			Bop(&q[0+i*n], temp, Bop_data);
			RNP::TBLAS::MultMV<'C'>(n, i+1, 1., z, n, temp, 1, 0., &rb[0+i*jmax], 1);
		}
		// The eigenvectors RA and RB  belonging to the found eigenvalues
		// are computed. The Schur vectors in VR and VS are replaced by the
		// eigenvectors of RA and RB
		RNP::GeneralizedEigensystem(k, ra, jmax, rb, jmax, alpha, beta, NULL, jmax, vsr, jmax, zwork, rwork);
		// Compute the eigenvectors belonging to the found eigenvalues
		// of A and put them in EIVEC
		RNP::TBLAS::MultMM<'N','N'>(n, k, k, 1., q, n, vsr, jmax, 0., eivec, ldeivec);
	}else{
		// Store the Schurvectors in eivec:
		//zcopy_(k * n, &work[q * n + 1], 1, &eivec[eivec_offset], 1);
		//RNP::TBLAS::Copy(k, aconv, 1, alpha, 1);
		//RNP::TBLAS::Copy(k, bconv, 1, beta, 1);
	}
	
	if(NULL == workspace){
		RNP::JDQZ_free(n, kmax, params, _work);
		delete _work;
	}
	return 0;
}


void RNP::JDQZ_alloc(size_t n, size_t k, const RNP::JDQZParams &params, RNP::JDQZWorkspace *workspace){
	if(NULL == workspace){ return; }
	size_t jmin = params.min_search_space;
	size_t jmax = params.max_search_space;
	jdqz_pick_search_space(k, jmin, jmax);
	workspace->zwork = new std::complex<double>[n*(2+5*jmax+3*k)+jmax*(10*jmax+5)];
	workspace->rwork = new double[8*jmax];
	workspace->iwork = new size_t[jmax];
}
void RNP::JDQZ_free(size_t n, size_t k, const RNP::JDQZParams &params, RNP::JDQZWorkspace *workspace){
	if(NULL == workspace){ return; }
	if(NULL != workspace->zwork){ delete [] workspace->zwork; workspace->zwork = NULL; }
	if(NULL != workspace->rwork){ delete [] workspace->rwork; workspace->rwork = NULL; }
	if(NULL != workspace->iwork){ delete [] workspace->iwork; workspace->iwork = NULL; }
}

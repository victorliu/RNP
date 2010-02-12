#ifndef _RNP_ITERATIVE_GENERALIZED_EIGENSYSTEMS_H_
#define _RNP_ITERATIVE_GENERALIZED_EIGENSYSTEMS_H_

#include "IterativeLinearSolvers.h"

namespace RNP{

struct JDQZWorkspace{
	std::complex<double> *zwork;
	double *rwork;
	size_t *iwork;
	JDQZWorkspace(){
		zwork = NULL;
		rwork = NULL;
		iwork = NULL;
	}
};

// Predefined eigenvalue sorters
bool JDQZSort_MaxReal(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
bool JDQZSort_MaxAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
bool JDQZSort_MinReal(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
bool JDQZSort_MinAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
bool JDQZSort_MaxImag(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
bool JDQZSort_MinImag(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
// JDQZSort_Nearest: data should be pointer to std::complex<double> that is the target
bool JDQZSort_Nearest(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);
// JDQZSort_NearestAbs: data should be pointer to double that is the target magnitude
bool JDQZSort_NearestAbs(const std::complex<double> &n1, const std::complex<double> &d1, const std::complex<double> &n2, const std::complex<double> &d2, void *data);

// Predefined linear solver adapters
struct JDQZ_LinearSolver_adapter_data{
	RNP::LinearSolverWorkspace *workspace;
	RNP::LinearSolverParams ls_params;
	JDQZ_LinearSolver_adapter_data(){
		workspace = NULL;
	}
};
void JDQZ_IDRs_adapter(
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
);
void JDQZ_GMRES_adapter(
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
);

// Predefined identity matrix adapter (for when B = I).
// Use this when possible since internal optimizations are used when
// it is known that B = I. The data parameter should point to size_t n.
void JDQZ_Identity(const std::complex<double> *, std::complex<double> *, void*);

// Important recommended settings:
//    Use testspace = 3 for targeted solves (1 or 2 also work)
//    Use testspace = 5 for everything else
//    testspace = 4 is not supported
struct JDQZParams{
	double eps;
	double lock;
	size_t min_search_space;
	size_t max_search_space;
	size_t max_iters;
	size_t max_mult;
	int testspace;
	// Testspace:
	// 1 = standard Petrov
	// 2 = Standard "variable" Petrov
	// 3 = harmonic Petrov
	// 4 = same as searchspace
	// 5 = B*v
	// 6 = A*v
	void (*linear_solver)(
		size_t n,
		// These next two arguments should be passed unmodified as the linear
		// solver's matrix-vector operator data. Failure to do so results in
		// undefined behavior.
		void (*Op)(const std::complex<double> *, std::complex<double> *, void*),
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
	void *linear_solver_data;
	JDQZParams(){
		eps = 1e-9;
		max_iters = 1000;
		lock = 1e-9;
		max_mult = 10000;
		// auto select
		min_search_space = 0;
		max_search_space = 0;
		testspace = 2;
	}
};

// Inputs:
//   n - size of the problem
//   k - number of eigenstuffs to compute
//   Aop - A routine that applies the A matrix to the first argument and stores it in the second
//   Bop - just like Aop
//   EvalSorter - sorter for eigenvalues (args: alpha1,beta1,alpha2,beta2,data)
//   preconditioner - Applies (approximation of) inv(A-tau*B) to its argument
// Ouptuts:
//   alpha,beta - numerator/denominator pairs for eigenvalues (length k)
//   vr,ldvr - eigenvalue matrix and leading dimension
//
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
	// User data parameters passed onto the corresponding arguments
	void *Aop_data = NULL,
	void *Bop_data = NULL,
	void *EvalSorter_data = NULL,
	void *preconditioner_data = NULL
);

void JDQZ_alloc(size_t n, size_t k, const JDQZParams &params, JDQZWorkspace *workspace);
void JDQZ_free(size_t n, size_t k, const JDQZParams &params, JDQZWorkspace *workspace);

}; // namespace RNP

#endif // _RNP_ITERATIVE_GENERALIZED_EIGENSYSTEMS_H_

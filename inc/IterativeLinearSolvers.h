#ifndef _RNP_ITERATIVE_LINEAR_SOLVERS_H_
#define _RNP_ITERATIVE_LINEAR_SOLVERS_H_

#include <complex>

// All solvers conform to this interface:
//
// int NAME(
// 	const size_t n,
// 	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
// 	std::complex<double>* b,
// 	std::complex<double>* x,
// 	LinearSolverParams &ls_params,
// 	const IDRs_params &params,
// 	LinearSolverWorkspace *workspace = NULL,
// 	void *Adata = NULL,
// 	// P is a precondition which simply solves P*x' = x,
// 	// where x i the first argument. The second parameter is user data,
// 	// which is passed in through Pdata.
// 	void (*P)(std::complex<double>*, void*) = NULL,
// 	void *Pdata = NULL
// );
//
// Arguments:
//
//  n          Dimension of the problem
//
//  A          A function which multiplies the matrix by the first argument
//             and returns the result in the second. The second argument must
//             be manually cleared. The third parameter is user data, passed
//             in through Adata.
//
//  b          RHS of the problem (assumed to be length n). Possibly
//             overwritten on exit.
//
//  x          The solution vector is returned here (assumed to be length n).
//
//  ls_params  Parameters common to all linear solvers. On exit, each field
//             is overwritten with the actual values of the solve.
//
//  params     Parameters specific to each solver.
//
//  workspace  Optional workspace.
//
//  Adata      See description of A.
//
//  P          A precondition which simply solves P*y = x, overwriting x with
//             y, where x is the first argument. The second parameter is user
//             data, which is passed in through Pdata.
//
//  Pdata      See description of P.
//
// Workspace: If not specified, the workspaces will be automatically
//            allocated internally. For best performance in repeated calls,
//            Use the NAME_alloc/NAME_free methods to preallocate the space.


namespace RNP{

struct LinearSolverWorkspace{
	std::complex<double> *zwork;
	double *rwork;
	size_t *iwork;
	LinearSolverWorkspace(){
		zwork = NULL;
		rwork = NULL;
		iwork = NULL;
	}
};

struct LinearSolverParams{
	size_t max_iterations; // max number of iterations; 0 = autoselect
	double tol;            // termination tolerance
	size_t max_mult;       // max number of matrix-vector multiplies allowed; 0 is unlimited
	bool x_initialized;    // if x contains an initial guess
	LinearSolverParams(){
		max_iterations = 0;
		tol = 1e-8;
		max_mult = 0;
		x_initialized = false;
	}
};

struct IDRs_params{
	size_t s; // the s in IDR(s); larger means more space needed
	double angle;
	IDRs_params(){
		s = 4;
		angle = 0.7;
	}
};


int IDRs(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	LinearSolverParams &ls_params,
	const IDRs_params &params,
	// Optional parameters
	LinearSolverWorkspace *workspace = NULL,
	//std::complex<double> *work = NULL, // if not null, should be length (3*n + 2*s)*(s+1)
	//size_t *iwork = NULL, // if not null, should be length params.s
	void *Adata = NULL,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*) = NULL,
	void *Pdata = NULL
);
void IDRs_alloc(size_t n, const IDRs_params &params, LinearSolverWorkspace *workspace);
void IDRs_free(size_t n, const IDRs_params &params, LinearSolverWorkspace *workspace);

struct GMRES_params{
	size_t m; // maximum degree of minimum residual polynomial
	GMRES_params(){
		m = 100;
	}
};

int GMRES(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	LinearSolverParams &ls_params,
	const GMRES_params &params,
	// Optional parameters
	LinearSolverWorkspace *workspace = NULL,
	//std::complex<double> *work = NULL, // if not null, should be length (3*n + 2*s)*(s+1)
	//double *iwork = NULL, // if not null, should be length params.s
	void *Adata = NULL,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*) = NULL,
	void *Pdata = NULL
);
void GMRES_alloc(size_t n, const GMRES_params &params, LinearSolverWorkspace *workspace);
void GMRES_free(size_t n, const GMRES_params &params, LinearSolverWorkspace *workspace);



struct BiCGSTABl_params{
	size_t l;
	double delta;
	double angle;
	BiCGSTABl_params(){
		l = 4;
		delta = 1e-2;
		angle = 0.7;
	}
};

int BiCGSTABl(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	LinearSolverParams &ls_params,
	const BiCGSTABl_params &params,
	// Optional parameters
	LinearSolverWorkspace *workspace = NULL,
	//std::complex<double> *work = NULL, // if not null, should be length (3*n + 2*s)*(s+1)
	//double *iwork = NULL, // if not null, should be length params.s
	void *Adata = NULL,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*) = NULL,
	void *Pdata = NULL
);
void BiCGSTABl_alloc(size_t n, const BiCGSTABl_params &params, LinearSolverWorkspace *workspace);
void BiCGSTABl_free(size_t n, const BiCGSTABl_params &params, LinearSolverWorkspace *workspace);

}; // namespace RNP

#endif // _RNP_ITERATIVE_LINEAR_SOLVERS_H_

#ifndef _RNP_IRA_H_
#define _RNP_IRA_H_

#include <complex>

// Implicitly Restarted Arnoldi method (ARPACK)
// BSD License

namespace RNP{
namespace IRA{
	typedef void (*Areg )(size_t n, const std::complex<double> *x, std::complex<double> *y, void *data);
	typedef void (*Bfunc)(size_t n, const std::complex<double> *x, std::complex<double> *y, void *data);
	
	// If a is more desireable than b, then SortFunc(a,b) should return false.
	typedef bool (*SortFunc)(const std::complex<double> &a, const std::complex<double> &b, void *data);
	// Predefined sort functions (data is not used)
	bool LargestMagnitude (const std::complex<double> &a, const std::complex<double> &b, void *data);
	bool SmallestMagnitude(const std::complex<double> &a, const std::complex<double> &b, void *data);
	bool LargestRealPart (const std::complex<double> &a, const std::complex<double> &b, void *data);
	bool SmallestRealPart(const std::complex<double> &a, const std::complex<double> &b, void *data);
	bool LargestImaginaryPart (const std::complex<double> &a, const std::complex<double> &b, void *data);
	bool SmallestImaginaryPart(const std::complex<double> &a, const std::complex<double> &b, void *data);
	
	struct Params{
		size_t max_iterations; // actual number of iterations is returned here
		double tol;            // actual tolerance used is returned here
		bool no_eigenvectors;  // if true, then v contains invariant subspace instead of eigenvectors
		Params(){
			max_iterations = 1000;
			tol = 0;                 // default to machine precision
			no_eigenvectors = false;
		}
	};
	struct Workspace{
		std::complex<double> *workd;  // length 3*n
		std::complex<double> *workl;  // length (3*n_arnoldi + 5)*n_arnoldi.
		std::complex<double> *workev; // length 2*n_arnoldi;
		std::complex<double> *resid;  // length n
		double *rwork;                // length n_arnoldi
		bool *bwork;                  // length n_arnoldi
	};
	
	bool AllocateWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace);
	void FreeWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace);
	
	// 1. If Bop is NULL, then we solve the non-generalized problem:
	//      A*x = lambda*x
	//    and Aop must perform: y <- A*x
	//    This is znaupd's mode 1.
	// 2. If Bop is not NULL, then we solve the generalized problem:
	//      A*x = lambda*B*x
	//    and Aop must perform: y <- inv(B)*A*x
	//    and Bop must perform: y <- B*x
	//    This is znaupd's mode 2.
	// The number of converged eigenpairs is returned.
	int Regular(
		size_t n, Areg Aop, Bfunc Bop,
		size_t n_wanted, size_t n_arnoldi,   // n_wanted+1 <= n_arnoldi <= n
		SortFunc sorter,
		std::complex<double> *w,             // eigenvalues (length n_wanted+1)
		std::complex<double> *v, size_t ldv, // eigenvectors (size ldv-by-n_arnoldi), ldv >= n
		Params *params = NULL,
		Workspace *workspace = NULL,
		void *Adata = NULL, void *Bdata = NULL, void *SortData = NULL);
	
	typedef void (*Ashifted)(size_t n, const std::complex<double> &shift, const std::complex<double> *x, std::complex<double> *y, void *data);
	// 1. If Bop is NULL, then we solve the non-generalized problem:
	//      A*x = lambda*x
	//    and Aop must perform: y <- inv(A - shift*I)*x
	//    This is znaupd's mode 3.
	// 2. If Bop is not NULL, then we solve the generalized problem:
	//      A*x = lambda*B*x
	//    and Aop must perform: y <- inv(A - shift*B)*x
	//    and Bop must perform: y <- B*x
	//    This is znaupd's mode 3.
	// The number of converged eigenpairs is returned.
	int ShiftInvert(
		size_t n, const std::complex<double> &shift, Ashifted Aop, Bfunc Bop,
		size_t n_wanted, size_t n_arnoldi,   // n_wanted+1 <= n_arnoldi <= n
		SortFunc sorter,
		std::complex<double> *w,             // eigenvalues (length n_wanted+1)
		std::complex<double> *v, size_t ldv, // eigenvectors (size ldv-by-n_arnoldi), ldv >= n
		Params *params = NULL,
		Workspace *workspace = NULL,
		void *Adata = NULL, void *Bdata = NULL, void *SortData = NULL);

}; // namespace IRA
}; // namespace RNP

#endif // _RNP_IRA_H_

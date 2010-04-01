#ifndef _RNP_LINEAR_SOLVE_H_
#define _RNP_LINEAR_SOLVE_H_

#include "TBLAS.h"
#include "TLASupport.h"

namespace RNP{

template <char trans='N'>
struct LinearSolve{
	template <class T>
	LinearSolve(size_t n, size_t nRHS, T *a, size_t lda, T *b, size_t ldb, int *info = NULL, size_t *pivots = NULL){
		if(NULL != info){ *info = 0; }
		if(0 == n || nRHS == 0){ return; }
		
		size_t *ipiv = pivots;
		if(NULL == pivots){
			ipiv = new size_t[n];
		}
		
		int ret = RNP::TLASupport::LUDecomposition(n, n, a, lda, ipiv);
		if(0 == ret){
			RNP::TLASupport::LUSolve<trans>(n, nRHS, a, lda, ipiv, b, ldb);
		}
		
		if(NULL == pivots){
			delete [] ipiv;
		}
		if(NULL != info){ *info = ret; }
	}
};

}; // namespace RNP

#ifdef RNP_HAVE_LAPACK
# include "LinearSolve_lapack.h"
#endif

#endif // _RNP_LINEAR_SOLVE_H_

#ifndef _RNP_LINEAR_SOLVE_LAPACK_H_
#define _RNP_LINEAR_SOLVE_LAPACK_H_

#include "LinearSolve.h"

extern "C" void zgesv_(const long int &n, const long int &nrhs, std::complex<double> *a, 
	const long int &lda, long int *ipiv, std::complex<double> *b, const long int &ldb, long int *info);

namespace RNP{

template <>
template <>
inline LinearSolve<'N'>::LinearSolve(size_t n, size_t nRHS, std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, int *info, size_t *pivots){
	if(NULL != info){ *info = 0; }
	if(0 == n || nRHS == 0){ return; }
	
	size_t *ipiv = pivots;
	if(NULL == pivots){
		ipiv = new size_t[n];
	}
	
	long int ret;
	zgesv_(n, nRHS, a, lda, (long int*)ipiv, b, ldb, &ret);
	
	if(NULL == pivots){
		delete [] ipiv;
	}
	if(NULL != info){ *info = ret; }
}

}; // namespace RNP

#endif // _RNP_LINEAR_SOLVE_LAPACK_H_

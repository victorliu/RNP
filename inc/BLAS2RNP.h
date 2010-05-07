#ifndef _BLAS_TO_RNP_H_
#define _BLAS_TO_RNP_H_

// Assumes that F2C was used to convert Fortran sources.
// Compile the C sources as C++.

#define _USE_MATH_DEFINES
#include <complex>
#include "TBLAS.h"

#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifdef __GNUC__
#define RNP_INLINE __attribute__((always_inline)) inline
#else
#define RNP_INLINE  inline
#endif

RNP_INLINE int RNP_FORTRAN_NAME(sswap,SSWAP)(long *n, float *x, long *incx, float *y, long *incy){
	RNP::TBLAS::Swap(*n, x, *incx, y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(dswap,DSWAP)(long *n, double *x, long *incx, double *y, long *incy){
	RNP::TBLAS::Swap(*n, x, *incx, y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(cswap,CSWAP)(long *n, singlecomplex *x, long *incx, singlecomplex *y, long *incy){
	RNP::TBLAS::Swap(*n, (std::complex<float>*)x, *incx, (std::complex<float>*)y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(zswap,ZSWAP)(long *n, doublecomplex *x, long *incx, doublecomplex *y, long *incy){
	RNP::TBLAS::Swap(*n, (std::complex<double>*)x, *incx, (std::complex<double>*)y, *incy); return 0;
}


RNP_INLINE int RNP_FORTRAN_NAME(sscal,SSCAL)(long *n, float &a, float *x, long *incx){
	RNP::TBLAS::Scale(*n, a, x, *incx); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(dscal,DSCAL)(long *n, double &a, double *x, long *incx){
	RNP::TBLAS::Scale(*n, a, x, *incx); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(cscal,CSCAL)(long *n, singlecomplex *a, singlecomplex *x, long *incx){
	RNP::TBLAS::Scale(*n, *(std::complex<float>*)a, (std::complex<float>*)x, *incx); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(zscal,ZSCAL)(long *n, doublecomplex *a, doublecomplex *x, long *incx){
	RNP::TBLAS::Scale(*n, *(std::complex<double>*)a, (std::complex<double>*)x, *incx); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(csscal,CSSCAL)(long *n, float *a, singlecomplex *x, long *incx){
	RNP::TBLAS::Scale(*n, *a, (std::complex<float>*)x, *incx); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(zdscal,ZDSCAL)(long *n, double *a, doublecomplex *x, long *incx){
	RNP::TBLAS::Scale(*n, *a, (std::complex<double>*)x, *incx); return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(scopy,SCOPY)(long *n, float *x, long *incx, float *y, long *incy){
	RNP::TBLAS::Copy(*n, x, *incx, y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(dcopy,DCOPY)(long *n, double *x, long *incx, double *y, long *incy){
	RNP::TBLAS::Copy(*n, x, *incx, y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(ccopy,CCOPY)(long *n, singlecomplex *x, long *incx, singlecomplex *y, long *incy){
	RNP::TBLAS::Copy(*n, (std::complex<float>*)x, *incx, (std::complex<float>*)y, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(zcopy,ZCOPY)(long *n, doublecomplex *x, long *incx, doublecomplex *y, long *incy){
	RNP::TBLAS::Copy(*n, (std::complex<double>*)x, *incx, (std::complex<double>*)y, *incy); return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(saxpy,SAXPY)(long *n, float *sa, float *sx, long *incx, float *sy, long *incy){
	RNP::TBLAS::Axpy(*n, *sa, sx, *incx, sy, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(daxpy,dAXPY)(long *n, double *sa, double *sx, long *incx, double *sy, long *incy){
	RNP::TBLAS::Axpy(*n, *sa, sx, *incx, sy, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(caxpy,CAXPY)(long *n, singlecomplex *sa, singlecomplex *sx, long *incx, singlecomplex *sy, long *incy){
	RNP::TBLAS::Axpy(*n, *(std::complex<float>*)sa, (std::complex<float>*)sx, *incx, (std::complex<float>*)sy, *incy); return 0;
}
RNP_INLINE int RNP_FORTRAN_NAME(zaxpy,ZAXPY)(long *n, doublecomplex *sa, doublecomplex *sx, long *incx, doublecomplex *sy, long *incy){
	RNP::TBLAS::Axpy(*n, *(std::complex<double>*)sa, (std::complex<double>*)sx, *incx, (std::complex<double>*)sy, *incy); return 0;
}

RNP_INLINE float RNP_FORTRAN_NAME(snrm2,SNRM2)(long *n, float *x, long *incx){
	return RNP::TBLAS::Norm2(*n, x, *incx);
}
RNP_INLINE double RNP_FORTRAN_NAME(dnrm2,DNRM2)(long *n, double *x, long *incx){
	return RNP::TBLAS::Norm2(*n, x, *incx);
}
RNP_INLINE float RNP_FORTRAN_NAME(scnrm2,SCNRM2)(long *n, singlecomplex *x, long *incx){
	return RNP::TBLAS::Norm2(*n, (std::complex<float>*)x, *incx);
}
RNP_INLINE double RNP_FORTRAN_NAME(dznrm2,DZNRM2)(long *n, doublecomplex *x, long *incx){
	return RNP::TBLAS::Norm2(*n, (std::complex<double>*)x, *incx);
}

RNP_INLINE float RNP_FORTRAN_NAME(sasum,SASUM)(long *n, const float *x, long *incx){
	return RNP::TBLAS::Asum(*n, x, *incx);
}
RNP_INLINE double RNP_FORTRAN_NAME(dasum,DASUM)(long *n, const double *x, long *incx){
	return RNP::TBLAS::Asum(*n, x, *incx);
}
RNP_INLINE float RNP_FORTRAN_NAME(scasum,SCASUM)(long *n, const singlecomplex *x, long *incx){
	return RNP::TBLAS::Asum(*n, (std::complex<float>*)x, *incx);
}
RNP_INLINE double RNP_FORTRAN_NAME(dzasum,DZASUM)(long *n, const doublecomplex *x, long *incx){
	return RNP::TBLAS::Asum(*n, (std::complex<double>*)x, *incx);
}

RNP_INLINE long int RNP_FORTRAN_NAME(isamax,ISAMAX)(long *n, const float *sx, long *incx){
	return RNP::TBLAS::MaximumIndex(*n, sx, *incx) + 1;
}
RNP_INLINE long int RNP_FORTRAN_NAME(idamax,IDAMAX)(long *n, const double *sx, long *incx){
	return RNP::TBLAS::MaximumIndex(*n, sx, *incx) + 1;
}
RNP_INLINE long int RNP_FORTRAN_NAME(icamax,ICAMAX)(long *n, const singlecomplex *sx, long *incx){
	return RNP::TBLAS::MaximumIndex(*n, (std::complex<float>*)sx, *incx) + 1;
}
RNP_INLINE long int RNP_FORTRAN_NAME(izamax,IZAMAX)(long *n, const doublecomplex *sx, long *incx){
	return RNP::TBLAS::MaximumIndex(*n, (std::complex<double>*)sx, *incx) + 1;
}


RNP_INLINE void RNP_FORTRAN_NAME(zdotc,ZDOTC)(doublecomplex *z, long *n, doublecomplex *x, long *incx, doublecomplex *y, long *incy){
	*(std::complex<double>*)z = RNP::TBLAS::ConjugateDot(*n, (std::complex<double>*)x, *incx, (std::complex<double>*)y, *incy);
}

RNP_INLINE int RNP_FORTRAN_NAME(zgeru,ZGERU)(long *m, long *n, doublecomplex *alpha, doublecomplex *x, long *incx, doublecomplex *y, long *incy, doublecomplex *a, long *lda){
	RNP::TBLAS::Rank1Update(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)x, *incx, (std::complex<double>*)y, *incy, (std::complex<double>*)a, *lda); return 0;
}


RNP_INLINE int RNP_FORTRAN_NAME(zgemv,ZGEMV)(char *trans, long *m, long *n, 
	doublecomplex *alpha, doublecomplex *a, long *lda, doublecomplex *x,
	long *incx, doublecomplex *beta, doublecomplex *y, long *incy, ftnlen trans_len = 0){
	if('N' == trans[0] || 'n' == trans[0]){
		RNP::TBLAS::MultMV<'N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)x, *incx, *(std::complex<double>*)beta, (std::complex<double>*)y, *incy);
	}else if('C' == trans[0] || 'c' == trans[0]){
		RNP::TBLAS::MultMV<'N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)x, *incx, *(std::complex<double>*)beta, (std::complex<double>*)y, *incy);
	}else if('T' == trans[0] || 't' == trans[0]){
		RNP::TBLAS::MultMV<'N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)x, *incx, *(std::complex<double>*)beta, (std::complex<double>*)y, *incy);
	}
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(zgemm,ZGEMM)(const char *transa, const char *transb, long *m, long *
	n, long *k, doublecomplex *alpha, doublecomplex *a, long *lda, 
	doublecomplex *b, long *ldb, doublecomplex *beta, doublecomplex *c, long *ldc, ftnlen transa_len = 0, ftnlen transb_len = 0){
	if('N' == transa[0] || 'n' == transa[0]){
		if('N' == transb[0] || 'n' == transb[0]){
			RNP::TBLAS::MultMM<'N','N'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('C' == transb[0] || 'c' == transb[0]){
			RNP::TBLAS::MultMM<'N','C'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('T' == transb[0] || 't' == transb[0]){
			RNP::TBLAS::MultMM<'N','T'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}
	}else if('C' == transa[0] || 'c' == transa[0]){
		if('N' == transb[0] || 'n' == transb[0]){
			RNP::TBLAS::MultMM<'C','N'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('C' == transb[0] || 'c' == transb[0]){
			RNP::TBLAS::MultMM<'C','C'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('T' == transb[0] || 't' == transb[0]){
			RNP::TBLAS::MultMM<'C','T'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}
	}else if('T' == transa[0] || 't' == transa[0]){
		if('N' == transb[0] || 'n' == transb[0]){
			RNP::TBLAS::MultMM<'T','N'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('C' == transb[0] || 'c' == transb[0]){
			RNP::TBLAS::MultMM<'T','C'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}else if('T' == transb[0] || 't' == transb[0]){
			RNP::TBLAS::MultMM<'T','T'>(*m, *n, *k, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb, *(std::complex<double>*)beta, (std::complex<double>*)b, *ldb);
		}
	}
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(ztrmm,ZTRMM)(char *side, char *uplo, char *transa, char *diag, long *m, long *n, doublecomplex *alpha, doublecomplex *a, long *lda, doublecomplex *b, long *ldb, ftnlen l1 = 0, ftnlen l2 = 0, ftnlen l3 = 0, ftnlen l4 = 0){
	if('L' == side[0] || 'l' == side[0]){
		if('L' == uplo[0] || 'l' == uplo[0]){
			if('U' == diag[0] || 'u' == diag[0]){
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'L','L','N','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'L','L','C','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'L','L','T','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}else{
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'L','L','N','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'L','L','C','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'L','L','T','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}
		}else{
			if('U' == diag[0] || 'u' == diag[0]){
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'L','U','N','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'L','U','C','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'L','U','T','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}else{
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'L','U','N','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'L','U','C','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'L','U','T','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}
		}
	}else{
		if('L' == uplo[0] || 'l' == uplo[0]){
			if('U' == diag[0] || 'u' == diag[0]){
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'R','L','N','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'R','L','C','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'R','L','T','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}else{
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'R','L','N','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'R','L','C','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'R','L','T','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}
		}else{
			if('U' == diag[0] || 'u' == diag[0]){
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'R','U','N','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'R','U','C','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'R','U','T','U'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}else{
				if('N' == transa[0] || 'n' == transa[0]){
					RNP::TBLAS::MultTrM<'R','U','N','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else if('C' == transa[0] || 'c' == transa[0]){
					RNP::TBLAS::MultTrM<'R','U','C','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}else{
					RNP::TBLAS::MultTrM<'R','U','T','N'>(*m, *n, *(std::complex<double>*)alpha, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
				}
			}
		}
	}
	return 0;
}

#endif // _BLAS_TO_RNP_H_

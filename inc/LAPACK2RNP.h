#ifndef _LAPACK_TO_RNP_H_
#define _LAPACK_TO_RNP_H_

// Assumes that F2C was used to convert Fortran sources.
// Compile the C sources as C++.

#define _USE_MATH_DEFINES
#include <complex>
#include "TLASupport.h"
#include <limits>

#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifdef __GNUC__
#define RNP_INLINE __attribute__((always_inline)) inline
#else
#define RNP_INLINE  inline
#endif

RNP_INLINE double RNP_FORTRAN_NAME(dlamch,DLAMCH)(char *s, ftnlen ls = 0){
	switch(s[0]){
	case 'S':
	case 's':
		return std::numeric_limits<double>::min();
	case 'P':
	case 'p':
		return 2*std::numeric_limits<double>::epsilon();
	case 'E':
	case 'e':
		return std::numeric_limits<double>::epsilon();
	case 'B':
	case 'b':
		return 2.0;
	case 'N':
	case 'n':
		return std::numeric_limits<double>::digits;
	case 'R':
	case 'r':
		return 1.0; // assume rounding occurs
	case 'M':
	case 'm':
		return std::numeric_limits<double>::min_exponent;
	case 'U':
	case 'u':
		return std::numeric_limits<double>::min(); // not quite right
	case 'L':
	case 'l':
		return std::numeric_limits<double>::max_exponent;
	case 'O':
	case 'o':
		return std::numeric_limits<double>::max();
	default:
		return 0;
	}
}

RNP_INLINE int RNP_FORTRAN_NAME(dlabad,DLABAD)(double *a, double *b){
	return 0;
}

RNP_INLINE double RNP_FORTRAN_NAME(dlapy2,dlapy2)(double *a, double *b){
	double aa = fabs(*a), ab = fabs(*b), r;
	if(aa > ab){
		r = ab/aa;
		return aa * sqrt(1. + r*r);
	}else{
		r = aa/ab;
		return ab * sqrt(1. + r*r);
	}
}

RNP_INLINE int RNP_FORTRAN_NAME(zlacpy,ZLACPY)(char *uplo, long *m, long *n, doublecomplex *a, long *lda, doublecomplex *b, long *ldb, ftnlen uplo_len = 0){
	if('L' == uplo[0] || 'l' == uplo[0]){
		RNP::TBLAS::CopyMatrix<'L'>(*m, *n, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
	}else if('U' == uplo[0] || 'u' == uplo[0]){
		RNP::TBLAS::CopyMatrix<'U'>(*m, *n, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
	}else{
		RNP::TBLAS::CopyMatrix<'A'>(*m, *n, (std::complex<double>*)a, *lda, (std::complex<double>*)b, *ldb);
	}
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(zlaset,ZLASET)(char *uplo, long *m, long *n, doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, long *lda, ftnlen uplo_len = 0){
	if('L' == uplo[0] || 'l' == uplo[0]){
		RNP::TBLAS::SetMatrix<'L'>(*m, *n, *(std::complex<double>*)alpha, *(std::complex<double>*)beta, (std::complex<double>*)a, *lda);
	}else if('U' == uplo[0] || 'u' == uplo[0]){
		RNP::TBLAS::SetMatrix<'L'>(*m, *n, *(std::complex<double>*)alpha, *(std::complex<double>*)beta, (std::complex<double>*)a, *lda);
	}else{
		RNP::TBLAS::SetMatrix<'L'>(*m, *n, *(std::complex<double>*)alpha, *(std::complex<double>*)beta, (std::complex<double>*)a, *lda);
	}
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(zlascl,ZLASCL)(char *mtype, long *kl, long *ku, double *cfrom, double *cto, long *m, long *n, doublecomplex *a, long *lda, long *info, ftnlen mtype_len = 0){
	switch(mtype[0]){
	case 'G':
	case 'g':
		RNP::TLASupport::RescaleMatrix<'G'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'L':
	case 'l':
		RNP::TLASupport::RescaleMatrix<'L'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'U':
	case 'u':
		RNP::TLASupport::RescaleMatrix<'U'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'H':
	case 'h':
		RNP::TLASupport::RescaleMatrix<'H'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'B':
	case 'b':
		RNP::TLASupport::RescaleMatrix<'B'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'Q':
	case 'q':
		RNP::TLASupport::RescaleMatrix<'Q'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	case 'Z':
	case 'z':
		RNP::TLASupport::RescaleMatrix<'Z'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	default:
		RNP::TLASupport::RescaleMatrix<'G'>(*kl, *ku, *cfrom, *cto, *m, *n, (std::complex<double>*)a, *lda); break;
	}
	return *info = 0;
}


RNP_INLINE int RNP_FORTRAN_NAME(zlarnv,ZLARNV)(long *idist, long *iseed, long *n, doublecomplex *x){
	int iiseed[4] = {iseed[0], iseed[1], iseed[2], iseed[3]};
	switch(*idist){
	case 1:
		RNP::TBLAS::RandomVector<1>(*n, (std::complex<double>*)x, iiseed); return 0;
	case 2:
		RNP::TBLAS::RandomVector<2>(*n, (std::complex<double>*)x, iiseed); return 0;
	case 3:
		RNP::TBLAS::RandomVector<3>(*n, (std::complex<double>*)x, iiseed); return 0;
	case 4:
		RNP::TBLAS::RandomVector<4>(*n, (std::complex<double>*)x, iiseed); return 0;
	case 5:
		RNP::TBLAS::RandomVector<5>(*n, (std::complex<double>*)x, iiseed); return 0;
	default:
		return 0;
	}
}

RNP_INLINE double RNP_FORTRAN_NAME(zlanhs,ZLANHS)(char *norm, long *n, doublecomplex *a, long *lda, double *work, ftnlen norm_len = 0){
	double ret = 0;
	switch(norm[0]){
	case 'M':
	case 'm':
		RNP::TLASupport::CheapHessenbergNorm<'M'>(*n, (std::complex<double>*)a, *lda, &ret, work);
	case '1':
	case 'O':
	case 'o':
		RNP::TLASupport::CheapHessenbergNorm<'1'>(*n, (std::complex<double>*)a, *lda, &ret, work);
	case 'I':
	case 'i':
		RNP::TLASupport::CheapHessenbergNorm<'I'>(*n, (std::complex<double>*)a, *lda, &ret, work);
	case 'F':
	case 'f':
	case 'E':
	case 'e':
		RNP::TLASupport::CheapHessenbergNorm<'F'>(*n, (std::complex<double>*)a, *lda, &ret, work);
	default:
		break;
	}
	return ret;
}

RNP_INLINE int RNP_FORTRAN_NAME(zlartg,ZLARTG)(doublecomplex *f, doublecomplex *g, double *cs, doublecomplex *sn, doublecomplex *r){
	RNP::TLASupport::GeneratePlaneRotation(*(std::complex<double>*)f, *(std::complex<double>*)g, cs, (std::complex<double>*)sn, (std::complex<double>*)r); return 0;
}

extern int ztrevc_(char howmny, bool *_select, size_t n, std::complex<double> *_t, size_t ldt, std::complex<double> *_vl, size_t ldvl, std::complex<double> *_vr, size_t ldvr, size_t mm, size_t *m, std::complex<double> *_work, double *rwork);
RNP_INLINE int RNP_FORTRAN_NAME(ztrevc,ZTREVC)(char *side, char *howmny, logical *select, long *n, doublecomplex *t, long *ldt, doublecomplex *vl, long *ldvl, doublecomplex *vr, long *ldvr, long *mm, long *m, doublecomplex *work, double *rwork, long *info, ftnlen l1 = 0, ftnlen l2 = 0){
	size_t sm;
	*info = ztrevc_(howmny[0], (bool *)select, *n, (std::complex<double>*)t, *ldt, (std::complex<double>*)vl, *ldvl, (std::complex<double>*)vr, *ldvr, *mm, &sm, (std::complex<double>*)work, rwork);
	*m = sm;
	return 0;
}

extern size_t zlahqr_(bool wantt, bool wantz, size_t n, size_t ilo, size_t ihi, std::complex<double> *h, int ldh, std::complex<double> *w, size_t iloz, size_t ihiz, std::complex<double> *z, size_t ldz);
RNP_INLINE int RNP_FORTRAN_NAME(zlahqr,ZLAHQR)(logical *wantt, logical *wantz, long *n, long *ilo, long *ihi, doublecomplex *h, long *ldh, doublecomplex *w, long *iloz, long *ihiz, doublecomplex *z, long *ldz, long *info){
	*info = zlahqr_(*wantt != 0, *wantz != 0, *n, *ilo, *ihi, (std::complex<double>*)h, *ldh, (std::complex<double>*)w, *iloz, *ihiz, (std::complex<double>*)z, *ldz);
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(zunm2r,ZUNM2R)(char *side, char *trans, long *m, long *n, long *k, doublecomplex *a, long *lda, doublecomplex *tau, doublecomplex *c, long *ldc, doublecomplex *work, long *info, ftnlen l1 = 0, ftnlen l2 = 0){
	if('L' == side[0] || 'l' == side[0]){
		if('C' == trans[0] || 'c' == trans[0]){
			RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'L','C'>(*m, *n, *k, (std::complex<double>*)a, *lda, (std::complex<double>*)tau, (std::complex<double>*)c, *ldc, (std::complex<double>*)work);
		}else{
			RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'L','N'>(*m, *n, *k, (std::complex<double>*)a, *lda, (std::complex<double>*)tau, (std::complex<double>*)c, *ldc, (std::complex<double>*)work);
		}
	}else{
		if('C' == trans[0] || 'c' == trans[0]){
			RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'R','C'>(*m, *n, *k, (std::complex<double>*)a, *lda, (std::complex<double>*)tau, (std::complex<double>*)c, *ldc, (std::complex<double>*)work);
		}else{
			RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'R','N'>(*m, *n, *k, (std::complex<double>*)a, *lda, (std::complex<double>*)tau, (std::complex<double>*)c, *ldc, (std::complex<double>*)work);
		}
	}
	return 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(zgeqr2,ZGEQR2)(long *m, long *n, doublecomplex *a, long *lda, doublecomplex *tau, doublecomplex *work, long *info){
	RNP::TLASupport::QRFactorization(*m, *n, (std::complex<double>*)a, *lda, (std::complex<double>*)tau, (std::complex<double>*)work);
	return *info = 0;
}

extern void ztrexc_(int n, std::complex<double> *t, int ldt, std::complex<double> *q, int ldq, int ifst, int ilst);
RNP_INLINE int RNP_FORTRAN_NAME(ztrexc,ZTREXC)(char *compq, long *n, doublecomplex *t, long *ldt, doublecomplex *q, long *ldq, long *ifst, long *ilst, long *info, ftnlen l1 = 0){
	ztrexc_(*n, (std::complex<double>*)t, *ldt, (std::complex<double>*)q, *ldq, *ifst-1, *ilst-1);
	return *info = 0;
}

RNP_INLINE int RNP_FORTRAN_NAME(ztrsen,ZTRSEN)(char *job, char *compq, logical *select, long *n, doublecomplex *t, long *ldt, doublecomplex *q, long *ldq, doublecomplex *w, long *m, double *s, double *sep, doublecomplex *work, long *lwork, long *info, ftnlen l1 = 0, ftnlen l2 = 0){
	if('N' == job[0] || 'n' == job[0]){
		if(-1 == *lwork){ *(std::complex<double>*)work = 1; return 0; }
		if(*m == *n || *m == 0){ return 0; }
		long ks = 0;
		long ierr;
		for(long k = 1; k <= *n; ++k){
			if(select[k-1]){ // Swap the K-th eigenvalue to position KS.
				++ks;
				if(k != ks){
					RNP_FORTRAN_NAME(ztrexc,ZTREXC)(compq, n, t, ldt, q, ldq, &k, &ks, &ierr);
				}
			}
		}
		for(long k = 0; k < *n; ++k){
			w[k] = t[k+k*(*ldt)];
		}
	}else{
		// not implemented
	}
	return 0;
}


#endif // _LAPACK_TO_RNP_H_

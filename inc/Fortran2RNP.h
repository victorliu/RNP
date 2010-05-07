#ifndef _FORTRAN_TO_RNP_H_
#define _FORTRAN_TO_RNP_H_

// Assumes that F2C was used to convert Fortran sources.
// Compile the C sources as C++.

#define _USE_MATH_DEFINES
#include <complex>

#define RNP_FORTRAN_NAME(LCASE,UCASE) LCASE ## _
#ifdef __GNUC__
#define RNP_INLINE __attribute__((always_inline)) inline
#else
#define RNP_INLINE  inline
#endif

RNP_INLINE void z_div(doublecomplex *z, doublecomplex *a, doublecomplex *b){
	*(std::complex<double>*)z = (*(std::complex<double>*)a) / (*(std::complex<double>*)b);
}
RNP_INLINE double d_imag(doublecomplex *z){
	return z->i;
}
RNP_INLINE void d_cnjg(doublecomplex *zc,doublecomplex *z){
	zc->r = z->r; zc->i = -z->i;
}


RNP_INLINE int s_copy(char *a, char *b, ftnlen la, ftnlen lb){

	char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend) *a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; ) *--aend = *--b;
#endif
	else{
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend) *a++ = *b++;
#ifndef NO_OVERWRITE
		else{
			a += lb;
			while(b < bend) *--a = *--bend;
			a += lb;
		}
#endif
		while(a < aend) *a++ = ' ';
	}
	return 0;
}
RNP_INLINE long s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb){
	char *a, *aend, *b, *bend;
	a = (char *)a0;
	b = (char *)b0;
	aend = a + la;
	bend = b + lb;

	if(la <= lb){
		while(a < aend)
			if(*a != *b) return( *a - *b );
			else{ ++a; ++b; }

		while(b < bend)
			if(*b != ' ') return( ' ' - *b );
			else ++b;
	}else{
		while(b < bend)
			if(*a == *b){ ++a; ++b; }
			else return( *a - *b );
		while(a < aend)
			if(*a != ' ') return(*a - ' ');
			else ++a;
	}
	return 0;
}

#endif // _FORTRAN_TO_RNP_H_

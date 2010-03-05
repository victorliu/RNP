#include "LinearSolve.h"
#include "Random.h"
#include "IO.h"
#include <complex>

typedef std::complex<double> complex_t;

int main(){
	const size_t n = 10;
	const size_t lda = n;
	complex_t *A = new complex_t[n*lda];
	complex_t *Acopy = new complex_t[n*lda];
	complex_t *b = new complex_t[n];
	complex_t *x = new complex_t[n];
	complex_t *r = new complex_t[n];
	
	for(size_t j = 0; j < n; ++j){
		double ra[2];
		for(size_t i = 0; i < n; ++i){
			RNP::Random::RandomRealsUniform01(2, ra);
			Acopy[i+j*lda] = A[i+j*lda] = complex_t(2*ra[0]-1, 2*ra[1]-1);
		}
		RNP::Random::RandomRealsUniform01(2, ra);
		x[j] = b[j] = complex_t(2*ra[0]-1, 2*ra[1]-1);
	}
	
	RNP::LinearSolve<'N'>(n,1, A,lda, x,n);
	
	RNP::TBLAS::MultMV<'N'>(n,n, 1., Acopy,lda, x,1, 0., r, 1);
	RNP::TBLAS::Axpy(n,-1.,b,1,r,1);
	
	std::cout << "Max residual: " << r[RNP::TBLAS::MaximumIndex(n, r, 1)] << std::endl;
	
	delete [] A;
	delete [] Acopy;
	delete [] b;
	delete [] x;
	delete [] r;
	return EXIT_SUCCESS;
}



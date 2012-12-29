#include <cstdlib>
#include <complex>
#include <iostream>
#include "Sparse.h"
#include "Random.h"
#include "LeastSquares.h"

typedef std::complex<double> complex_t;
typedef RNP::Sparse::TCCSMatrix<complex_t> ccsmatrix;

void Aop(bool trans, std::complex<double> *x, std::complex<double> *y, void *data){ 
	const ccsmatrix *A = reinterpret_cast<ccsmatrix*>(data);
	if(trans){
		RNP::Sparse::MultMV<'C'>(*A, y, x, complex_t(1), complex_t(1));
	}else{
		RNP::Sparse::MultMV<'N'>(*A, x, y, complex_t(1), complex_t(1));
	}
}

int main(){
	srand(time(0));
	size_t m = 1000;
	size_t n = 500;
	std::complex<double> *b = new std::complex<double>[m];
	std::complex<double> *x = new std::complex<double>[n];
	std::complex<double> *r = new std::complex<double>[m];
	ccsmatrix::entry_map_t entries;
	for(size_t i = 0; i < n; ++i){
		entries[ccsmatrix::index_t(i,i)] = 1;
	}
	for(size_t i = 0; i < m; ++i){
		double ra[2];
		RNP::Random::RandomRealsUniform01(2, ra);
		b[i] = std::complex<double>(ra[0], 2*ra[1]-1);
	}
	
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < 4; ++j){
			int row = rand()%m;
			int col = i;
			double ra[2];
			RNP::Random::RandomRealsUniform01(2, ra);
			entries[ccsmatrix::index_t(row,col)] += std::complex<double>(2*ra[0]-1, 2*ra[1]-1);
		}
	}
	ccsmatrix A(m, n, entries);

	//print(&A); std::cout << std::endl << std::endl;
	//print(n, b); std::cout << std::endl;
	
	RNP::LSQR<complex_t>::Params params;
	RNP::LSQR<complex_t>::OptionalOutputs out;
	RNP::LSQR<complex_t>(m, n, &Aop, complex_t(0), b, x, params, NULL, (void*)&A, &out);
	
	std::cout << "info   = " << out.info  << std::endl;
	std::cout << "anorm  = " << out.anorm << std::endl;
	std::cout << "acond  = " << out.acond << std::endl;
	std::cout << "rnorm  = " << out.rnorm << std::endl;
	std::cout << "arnorm = " << out.arnorm << std::endl;
	std::cout << "xnorm  = " << out.xnorm << std::endl;
	
	delete [] b;
	delete [] x;
	delete [] r;
	return 0;
}

#include "Sparse.h"
#include "Random.h"
#include "IterativeLinearSolvers.h"
#include "TBLAS.h"
#include <complex>
#include <iostream>

typedef RNP::Sparse::TCCSMatrix<std::complex<double> > ccsmatrix;

void Aop(const std::complex<double> *x, std::complex<double> *y, void *data){ 
	const ccsmatrix *A = reinterpret_cast<ccsmatrix*>(data);
	RNP::Sparse::MultMV(*A, x, y);
}

int main(){
	srand(time(0));
	size_t n = 1000;
	std::complex<double> *b = new std::complex<double>[n];
	std::complex<double> *x = new std::complex<double>[n];
	std::complex<double> *r = new std::complex<double>[n];
	ccsmatrix::entry_map_t entries;
	for(size_t i = 0; i < n; ++i){
		entries[ccsmatrix::index_t(i,i)] = 1;
		double ra[2];
		RNP::Random::RandomRealsUniform01(2, ra);
		b[i] = std::complex<double>(ra[0], 2*ra[1]);
	}
	
	for(size_t i = 0; i < n; ++i){
		int row = rand()%n;
		int col = rand()%n;
		double ra[2];
		RNP::Random::RandomRealsUniform01(2, ra);
		entries[ccsmatrix::index_t(row,col)] += std::complex<double>(2*ra[0]-1, 2*ra[1]-1);
	}
	ccsmatrix A(n, entries);

	//print(&A); std::cout << std::endl << std::endl;
	//print(n, b); std::cout << std::endl;
	
	RNP::LinearSolverParams ls_params;
	RNP::IDRs_params params;
	RNP::IDRs(n, &Aop, b, x, ls_params, params, NULL, (void*)&A);
	
	//print(n, x); std::cout << std::endl;
	
	RNP::Sparse::MultMV(A, x, r);
	for(size_t i = 0; i < n; ++i){
		r[i] -= b[i];
	}
	double normr = RNP::TBLAS::Norm2(n, r, 1);
	double normb = RNP::TBLAS::Norm2(n, b, 1);
	std::cout << "normb = " << normb << std::endl;
	std::cout << "err = " << normr/normb << std::endl;
	
	delete [] b;
	delete [] x;
	delete [] r;
	return 0;
}

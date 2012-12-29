#include <cstdlib>
#include "Sparse.h"
#include "Random.h"
#include "IterativeLinearSolvers.h"
#include "TBLAS.h"
#include "IO.h"
#include <complex>
#include <iostream>

#include "IterativeGeneralizedEigensystems.h"

typedef RNP::Sparse::TCCSMatrix<std::complex<double> > ccsmatrix;

void Aop(const std::complex<double> *x, std::complex<double> *y, void *data){ 
	const ccsmatrix *A = reinterpret_cast<ccsmatrix*>(data);
	RNP::Sparse::MultMV<'N'>(*A, x, y);
}
void Bop(const std::complex<double> *x, std::complex<double> *y, void *data){ 
	const ccsmatrix *B = reinterpret_cast<ccsmatrix*>(data);
	RNP::Sparse::MultMV<'N'>(*B, x, y);
}
void Precond(std::complex<double> *x, void *data){
	size_t *n = reinterpret_cast<size_t*>(data);
}

void jdqz_manual_test(){
	const size_t n = 100;
	const size_t k = 5;
	std::complex<double> *alpha = new std::complex<double>[k];
	std::complex<double> *beta  = new std::complex<double>[k];
	std::complex<double> *resid = new std::complex<double>[n];
	std::complex<double> *temp = new std::complex<double>[n];
	std::complex<double> *eivec = new std::complex<double>[n*k];
	ccsmatrix::entry_map_t entriesA, entriesB;
	for(size_t i = 0; i < n; ++i){
		entriesA[ccsmatrix::index_t(i,i)] = i+1.0;
		entriesB[ccsmatrix::index_t(i,i)] = 1.0/(i+1.0);
	}
	ccsmatrix A(n, entriesA);
	ccsmatrix B(n, entriesB);
	
	std::complex<double> target(31.0);
	RNP::JDQZParams params;
	params.linear_solver = &RNP::JDQZ_GMRES_adapter;
	RNP::JDQZ_LinearSolver_adapter_data linsolve_data;
	params.linear_solver_data = &linsolve_data;
	params.testspace = 3;
	JDQZ(n, k, &Aop, &Bop, alpha, beta, eivec, n, target, &RNP::JDQZSort_Nearest, &Precond, params, NULL, (void*)&A, (void*)&B, (void*)&target, (void*)&n);
	
	for(size_t i = 0; i < k; ++i){
		Aop(&eivec[0+i*n], resid, (void*)&A);
		RNP::TBLAS::Scale(n, beta[i], resid, 1);
		Bop(&eivec[0+i*n], temp, (void*)&B);
		RNP::TBLAS::Axpy(n, -alpha[i], temp, 1, resid, 1);
		std::cout << "lambda[" << i << "] = ";
		RNP::IO::Print(alpha[i]/beta[i]);
		std::cout << std::endl;
		std::cout << "|beta*A*x - alpha*B*x| = ";
		RNP::IO::Print(RNP::TBLAS::Norm2(n, resid, 1));
		std::cout << std::endl;
	}
	delete [] alpha;
	delete [] beta;
	delete [] resid;
	delete [] temp;
	delete [] eivec;
}

int main(){
	//srand((unsigned int)time(0));
	srand(1);
	const size_t n = 100;
	const size_t k = 5;
	std::complex<double> *alpha = new std::complex<double>[k];
	std::complex<double> *beta  = new std::complex<double>[k];
	std::complex<double> *resid = new std::complex<double>[n];
	std::complex<double> *temp = new std::complex<double>[n];
	std::complex<double> *eivec = new std::complex<double>[n*k];
	ccsmatrix::entry_map_t entriesA, entriesB;
	for(size_t i = 0; i < n; ++i){
		entriesA[ccsmatrix::index_t(i,i)] = 1;
		entriesB[ccsmatrix::index_t(i,i)] = 1;
	}
	for(size_t i = 0; i < n; ++i){
		int row = rand()%n;
		int col = rand()%n;
		double ra[2];
		RNP::Random::RandomRealsUniform01(2, ra);
		//entriesA[ccsmatrix::index_t(row,col)] += std::complex<double>(2.*ra[0]-1., 2.*ra[1]-1.);
		RNP::Random::RandomRealsUniform01(2, ra);
		//entriesB[ccsmatrix::index_t(row,col)] += std::complex<double>(2.*ra[0]-1., 2.*ra[1]-1.);
	}
	ccsmatrix A(n, entriesA);
	ccsmatrix B(n, entriesB);

	//print(&A); std::cout << std::endl << std::endl;
	//print(n, b); std::cout << std::endl;
	
	std::complex<double> target(0);
	RNP::JDQZParams params;
	params.linear_solver = &RNP::JDQZ_GMRES_adapter;
	RNP::JDQZ_LinearSolver_adapter_data linsolve_data;
	params.linear_solver_data = &linsolve_data;
	params.testspace = 5;
	RNP::JDQZ(n, k, &Aop, &Bop, alpha, beta, eivec, n, target, &RNP::JDQZSort_MinReal, &Precond, params, NULL, (void*)&A, (void*)&B, NULL, (void*)&n);
	
	for(size_t i = 0; i < k; ++i){
		Aop(&eivec[0+i*n], resid, (void*)&A);
		RNP::TBLAS::Scale(n, beta[i], resid, 1);
		Bop(&eivec[0+i*n], temp, (void*)&B);
		RNP::TBLAS::Axpy(n, -alpha[i], temp, 1, resid, 1);
		std::cout << "lambda[" << i << "] = ";
		RNP::IO::Print(alpha[i]/beta[i]);
		std::cout << std::endl;
		std::cout << "|beta*A*x - alpha*B*x| = ";
		RNP::IO::Print(RNP::TBLAS::Norm2(n, resid, 1));
		std::cout << std::endl;
	}
	delete [] alpha;
	delete [] beta;
	delete [] resid;
	delete [] temp;
	delete [] eivec;
	return 0;
}

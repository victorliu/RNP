#include "RootFinding.h"
#include <cstdlib>
#include <iostream>
#include <ctime>

typedef RNP::RootFinding::Mueller::complex_t complex_t;

double frand(){
	return (double)rand()/(double)RAND_MAX - 0.5;
}

complex_t f(const complex_t &z, void *data){
	static const complex_t r[5] = {
		complex_t(0.3, 0.5),
		complex_t(0.2, 1.5),
		complex_t(-0.1, 0.9),
		complex_t(10, 0.01),
		complex_t(10.01, 0.02)
	};
	static const complex_t a[5] = {
		complex_t(2, 5),
		complex_t(1, 0.2),
		complex_t(-6, 1),
		complex_t(-2, 10),
		complex_t(1, 10)
	};
	complex_t ret(1.0);
	for(size_t i = 0; i < 5; ++i){
		ret += a[i]/(z - r[i]);
	}
	ret = 1./ret;
	std::cout << "f" << z << " = " << ret << std::endl;
	return ret;
}
#define NSEEK 3
int main(){
	srand(1);
	RNP::RootFinding::Mueller::Params params;
	RNP::RootFinding::Mueller::ComplexRoot root[NSEEK];
	params.tolerance = 0;
	params.max_evaluations = 0;
	for(size_t i = 0; i < NSEEK; ++i){
		root[i].z = complex_t(frand(), frand());
		root[i].dz = 0;
	}
	root[0].z = complex_t(0.29,0.51);
	RNP::RootFinding::Mueller::FindRoots(&f, NSEEK, root, &params, NULL);
	for(size_t i = 0; i < NSEEK; ++i){
		if(root[i].converged){
			std::cout << "Converged: ";
		}else{
			std::cout << "Unconverged: ";
		}
		std::cout << root[i].z << std::endl;
	}
	return 0;
}

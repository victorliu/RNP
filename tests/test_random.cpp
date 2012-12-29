#include "Random.h"
#include <iostream>
#include <cstdlib>

int main(){
	std::cout << "All number should appear random in (0,1)" << std::endl;
	std::cout << "Testing default seed, one at a time:" << std::endl;
	for(size_t i = 0; i < 100; ++i){
		double x;
		RNP::Random::RandomRealsUniform01(1, &x);
		std::cout << x << std::endl;
	}
	std::cout << "Testing default seed, 100 in a single go:" << std::endl;
	double *v = new double[100];
	RNP::Random::RandomRealsUniform01(100, v);
	for(size_t i = 0; i < 100; ++i){
		std::cout << v[i] << std::endl;
	}
	int seed[4] = {time(0)%4096, time(0)%4096, time(0)%4096, (time(0)%4096)|1};
	std::cout << "Testing time(0) seed, one at a time:" << std::endl;
	for(size_t i = 0; i < 100; ++i){
		double x;
		RNP::Random::RandomRealsUniform01(1, &x, seed);
		std::cout << x << std::endl;
	}
	std::cout << "Testing default seed, 100 in a single go:" << std::endl;
	RNP::Random::RandomRealsUniform01(100, v, seed);
	for(size_t i = 0; i < 100; ++i){
		std::cout << v[i] << std::endl;
	}
	delete [] v;
	return EXIT_SUCCESS;
}

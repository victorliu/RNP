#include <IO.h>

int main(){
	double v = 153.234;
	std::complex<double> c(45.8,-0.24);
	RNP::IO::Print(v) << std::endl;
	RNP::IO::Print(c) << std::endl;
	return 0;
}

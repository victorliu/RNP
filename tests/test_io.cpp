#include "IO.h"

int main(){
	double v = 153.234;
	std::complex<double> c(45.8,-0.24);
	Print(v) << std::endl;
	Print(c) << std::endl;
	return 0;
}

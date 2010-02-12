% RNP Documentation -- IO
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# RNP::IO

The IO module contains convenience functions for outputting numbers, vectors, and matrices to C++ STL streams.
The output format can be tailored for easy import into numerical packages.
Currently, there are preprocessor flags `RNP_OUTPUT_MATHEMATICA` and `RNP_OUTPUT_MATLAB` for the respective software programs.
All functions are in namespace `RNP::IO`.

Note that Mathematica output does not convert the exponential `E` in floating point numbers into `*^`; this conversion must be performed manually.

## Contents

* [Print](#IO_Print)
* [PrintVector](#IO_PrintVector)
* [PrintMatrix](#IO_PrintMatrix)

---
## Print<a name="IO_Print" />

Outputs a number.
Returns the output stream passed in.

### Prototype

	template <class T>
	std::ostream& Print(const T &value, std::ostream &os = std::cout);

### Arguments:

=value=
	The scalar value to output.
=os=
	The output stream to use. Defaults to `std::cout`.

---
## PrintVector<a name="IO_PrintVector" />

Outputs a vector.
Returns the output stream passed in.

### Prototype

	template <class T>
	std::ostream& PrintVector(size_t n, const T *x, size_t incx = 1, std::ostream &os = std::cout);

### Arguments:

=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.
=os=
	The output stream to use. Defaults to `std::cout`.




---
## PrintMatrix<a name="IO_PrintMatrix" />

Outputs a rectangular matrix.
Returns the output stream passed in.

### Prototype

	template <class T>
	std::ostream& PrintMatrix(size_t m, size_t n, const T *a, size_t lda, std::ostream &os = std::cout);

### Arguments:

=m=
	Number of rows of the matrix.
=n=
	Number of columns of the matrix.
=a=
	Pointer to the first element of the matrix.
=lda=
	Leading dimension of the matrix (typically the number of rows of `a` unless `a` is a submatrix of a larger matrix).
=os=
	The output stream to use. Defaults to `std::cout`.


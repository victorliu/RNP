% RNP Documentation -- TBLAS
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# RNP::TBLAS

TBLAS is a templated BLAS, with some design modifications.
All functions are in namespace `RNP::TBLAS`.
Optionally, if the RNP Random module is used, define `RNP_TBLAS_USE_RANDOM` to enable functions for filling in vectors or matrices with random values.

## Contents

### BLAS Level 1 functions
* [Fill](#TBLAS_Fill)
* [Swap](#TBLAS_Swap)
* [Scale](#TBLAS_Scale)
* [Copy](#TBLAS_Copy)
* [Axpy](#TBLAS_Axpy)
* [Dot](#TBLAS_Dot)
* [Norm2](#TBLAS_Norm2)
* [Asum](#TBLAS_Asum)
* [MaximumIndex](#TBLAS_MaximumIndex)

### BLAS Level 2 functions
* [Axpy](#TBLAS_Axpy_mat)
* [MultMV](#TBLAS_MultMV)
* [MultTrV](#TBLAS_MultTrV)
* [SolveTrV](#TBLAS_SolveTrV)
* [Rank1Update](#TBLAS_Rank1Update)
* [ConjugateRank1Update](#TBLAS_ConjugateRank1Update)

### BLAS Level 3 functions
* [MultMM](#TBLAS_MultMM)
* [SolveTrM](#TBLAS_SolveTrM)

### LAPACK utility functions
* [Conjugate](#TBLAS_Conjugate)
* [CopyMatrix](#TBLAS_CopyMatrix)
* [SetMatrix](#TBLAS_SetMatrix)
* [RandomVector](#TBLAS_RandomVector)


---
## Fill<a name="TBLAS_Fill" />

Fills a vector with a constant value.

This function does not correspond to any BLAS function.

### Prototype

	template <class TV, class T>
	void Fill(size_t n, const TV &value, T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=value=
	The value each element will be set to.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.

---
## Swap<a name="TBLAS_Swap" />

Swaps two vectors (`x <--> y`).

This function corresponds to BLAS functions `zswap`, `dswap`, `cswap`, and `sswap`.

### Prototype

	template <class T>
	void Swap(size_t n, T *x, size_t incx, T *y, size_t incy);

### Arguments:

=n=
	Length of the vectors.
=x=
	Pointer to the first element of vector `x`.
=incx=
	Increment (stride) of vector `x`.
=y=
	Pointer to the first element of vector `y`.
=incy=
	Increment (stride) of vector `y`.


---
## Scale<a name="TBLAS_Scale" />

Scales a vector `x` by `alpha` (`x <-- alpha * x`).

This function corresponds to BLAS functions `zscal`, `dscal`, `cscal`, and `sscal`.

### Prototype

	template <class TA, class T>
	void Scale(size_t n, const TA &alpha, T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=alpha=
	The scale factor.
=x=
	Pointer to the first element of vector `x`.
=incx=
	Increment (stride) of vector `x`.

---
## Copy<a name="TBLAS_Copy" />

Copies one vector into another (`dst <-- src`).

This function corresponds to BLAS functions `zcopy`, `dcopy`, `ccopy`, and `scopy`.

### Prototype

	template <class T>
	void Copy(size_t n, const T *src, size_t incsrc, T *dst, size_t incdst);

### Arguments:

=n=
	Length of the vectors.
=src=
	Pointer to the first element of the source vector.
=incsrc=
	Increment (stride) of source vector.
=dst=
	Pointer to the first element of the destination vector.
=incdst=
	Increment (stride) of destination vector.


---
## Axpy<a name="TBLAS_Axpy" />

Performs the scaled addition of one vector to another (`y <-- alpha * x + y`).

This function corresponds to BLAS functions `zaxpy`, `daxpy`, `caxpy`, and `saxpy`.

### Prototype

	template <class S, class T>
	void Axpy(size_t n, const S &alpha, const T *x, size_t incx, T *y, size_t incy);

### Arguments:

=n=
	Length of the vectors.
=alpha=
	The scale factor that is applied to source vector `x`.
=x=
	Pointer to the first element of the source vector `x`.
=incx=
	Increment (stride) of source vector `x`.
=y=
	Pointer to the first element of the destination vector `y`.
=incy=
	Increment (stride) of destination vector `y`.



---
## Dot<a name="TBLAS_Dot" />

Returns the dot product (scalar product) of two vectors (`x' * y`).

This function corresponds to BLAS functions `zdot`, `ddot`, `cdot`, and `sdot`.

### Prototype

	template <class T>
	T Dot(size_t n, const T *x, size_t incx, T *y, size_t incy);

### Arguments:

=n=
	Length of the vectors.
=x=
	Pointer to the first element of the vector `x`.
=incx=
	Increment (stride) of vector `x`.
=y=
	Pointer to the first element of the vector `y`.
=incy=
	Increment (stride) of vector `y`.





---
## ConjugateDot<a name="TBLAS_ConjugateDot" />

Returns the Hermitian dot product (scalar product) of two vectors (`x^H * y`).

This function corresponds to BLAS functions `zdotc` and `cdotc`.

### Prototype

	template <class T>
	T ConjugateDot(size_t n, const T *x, size_t incx, T *y, size_t incy);

### Arguments:

=n=
	Length of the vectors.
=x=
	Pointer to the first element of the vector `x`. This is the vector that is conjugated.
=incx=
	Increment (stride) of vector `x`.
=y=
	Pointer to the first element of the vector `y`.
=incy=
	Increment (stride) of vector `y`.


---
## Norm2<a name="TBLAS_Norm2" />

Returns the Euclidean norm (2-norm) of a vector (`sqrt(x^H * x)`).
The return type is the corresponding real-number type of template type `T` (i.e. if `T=double` then `real_type=double`, and if `T=std::complex<double>` then `real_type=double`).

This function corresponds to BLAS functions `dznrm2`, `dnrm2`, `scnrm2`, and `snrm2`.

### Prototype

	template <class T>
	real_type Norm2(size_t n, const T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.

---
## Asum<a name="TBLAS_Asum" />

Returns the sum of magnitudes of a vector (kind of like a 1-norm).
For complex numbers, the magnitude is taken to be the sum of the absolute values of the real and imaginary parts (the 1-norm of its vector representation in the complex plane). This is meant to be a cheap norm-like quantity.
The return type is the corresponding real-number type of template type `T` (i.e. if `T=double` then `real_type=double`, and if `T=std::complex<double>` then `real_type=double`).

This function corresponds to BLAS functions `dzasum`, `dasum`, `scasum`, and `sasum`.

### Prototype

	template <class T>
	typename _RealOrComplexChooser<T>::real_type Asum(size_t n, const T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.



---
## MaximumIndex<a name="TBLAS_MaximumIndex" />

Returns the (0-based) index of the element of a vector with the largest magnitude.
For complex numbers, the magnitude is taken to be the sum of the absolute values of the real and imaginary parts (the 1-norm of its vector representation in the complex plane).

This function corresponds to BLAS functions `izamax`, `idamax`, `icamax`, and `isamax`.

### Prototype

	template <class T>
	size_t MaximumIndex(size_t n, const T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.




---
## Axpy<a name="TBLAS_Axpy_mat" />

Performs the scaled addition of one matrix to another (`y <-- alpha * x + y`).

This function does not correspond to any BLAS function.

### Prototype

	template <class S, class T>
	void Axpy(size_t m, size_t n, const S &alpha, const T *x, size_t ldx, T *y, size_t ldy);

### Arguments:

=m=
	Number of rows of the matrices `x` and `y`.
=n=
	Number of columns of the matrices.
=alpha=
	The scale factor applied to `x`.
=x=
	Pointer to the first element of matrix `x`.
=ldx=
	Leading dimension of matrix `x` (typically the number of rows of `x` unless `x` is a submatrix of a larger matrix).
=y=
	Pointer to the first element of vector `y`.
=incy=
	Increment (stride) of the vector `y`.



---
## MultMV<a name="TBLAS_MultMV" />

Performs a matrix-vector multiply (`y <-- alpha * op(A) * x + beta * y`, where `op` is specified by template parameter `trans`).
No index or dimension checking is performed.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameter corresponds to the classical `BLAS` parameter `trans` and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::MultMV<'N'>(m, n, alpha, a, lda, x, incx, beta, y, incy);

This function corresponds to BLAS functions `zgemv`, `dgemv`, `cgemv`, and `sgemv`.

### Prototype

	template <char trans>
	struct MultMV{
		template <class A, class B, class T>
		MultMV(size_t m, size_t n, const A &alpha, const T *a, size_t lda,
		       const T *x, size_t incx,
		       const B &beta, T *y, size_t incy);
	};

### Arguments:

=trans=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `trans=N`, `op(A) = A`.
	* If `trans=T`, `op(A) = A'` (transpose of A)
	* If `trans=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=m=
	Number of rows of `op(A)` (i.e. if `trans=T`, `m` is the number of columns of `A`).
=n=
	Number of columns of `op(A)`.
=alpha=
	The scale factor applied to `op(A) * x` (for a simple matrix-vector multiply, set this to 1).
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=x=
	Pointer to the first element of vector `x`.
=incx=
	Increment (stride) of the vector `x`.
=beta=
	The scale factor applied to `y` before adding the product (for a simple matrix-vector multiply, set this to 0).
=y=
	Pointer to the first element of destination vector `y`.
=incy=
	Increment (stride) of the destination vector `y`.


---
## MultTrV<a name="TBLAS_MultTrV" />

Performs a matrix-vector multiply (`x <-- op(A) * x`, where `op` is specified by template parameter `trans`) where the matrix `A` is a triangular portion of a rectangular matrix.
No index or dimension checking is performed.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameters correspond to the classical `BLAS` parameters `uplo`, `trans`, and `diag`, and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::MultTrV<'L','N','N'>(n, a, lda, x, incx);

This function corresponds to BLAS functions `ztrmv`, `dtrmv`, `ctrmv`, and `strmv`.

### Prototype

	template <char uplo, char trans, char diag>
	struct MultTrV{
		template <class T>
		MultTrV(size_t n, const T *a, size_t lda, T *x, size_t incx);
	};

### Arguments:

=uplo=
	Specifies whether `A` is upper or lower triangular.
	This is a template parameter which must be explicitly specified.
	* If `uplo=U`, `A` is an upper triangular matrix.
	* If `uplo=L`, `A` is a lower triangular matrix.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=trans=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `trans=N`, `op(A) = A`.
	* If `trans=T`, `op(A) = A'` (transpose of A)
	* If `trans=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=diag=
	Whether the diagonal of `A` is all ones (unit triangular).
	This is a template parameter which must be explicitly specified.
	* If `diag=U`, `A` is assumed to be unit triangular.
	* If `diag=N`, `A` is not assumed to be unit triangular.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=n=
	Number of rows/columns of `op(A)`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=x=
	Pointer to the first element of vector `x`.
=incx=
	Increment (stride) of the vector `x`.


---
## SolveTrV<a name="TBLAS_SolveTrV" />

Performs a triangular matrix solve (`x <-- inv(op(A)) * x`, where `op` is specified by template parameter `trans`) where the matrix `A` is a triangular portion of a rectangular matrix.
No index or dimension checking is performed.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameters correspond to the classical `BLAS` parameters `uplo`, `trans`, and `diag`, and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::SolveTrV<'L','N','N'>(n, a, lda, x, incx);

This function corresponds to BLAS functions `ztrsv`, `dtrsv`, `ctrsv`, and `strsv`.

### Prototype

	template <char uplo, char trans, char diag>
	struct SolveTrV{
		template <class T>
		SolveTrV(size_t n, const T *a, size_t lda, T *x, size_t incx);
	};

### Arguments:

=uplo=
	Specifies whether `A` is upper or lower triangular.
	This is a template parameter which must be explicitly specified.
	* If `uplo=U`, `A` is an upper triangular matrix.
	* If `uplo=L`, `A` is a lower triangular matrix.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=trans=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `trans=N`, `op(A) = A`.
	* If `trans=T`, `op(A) = A'` (transpose of A)
	* If `trans=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=diag=
	Whether the diagonal of `A` is all ones (unit triangular).
	This is a template parameter which must be explicitly specified.
	* If `diag=U`, `A` is assumed to be unit triangular.
	* If `diag=N`, `A` is not assumed to be unit triangular.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=n=
	Number of rows/columns of `op(A)`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=x=
	Pointer to the first element of vector `x`.
=incx=
	Increment (stride) of the vector `x`.


---
## Rank1Update<a name="TBLAS_Rank1Update" />

Performs a rank 1 update of a matrix (`A <-- alpha * x * y' + A`.
No index or dimension checking is performed.

This function corresponds to BLAS functions `zgeru`, `dgeru`, `cgeru`, and `sgeru`.

### Prototype

	template <class A, class T>
	void Rank1Update(size_t m, size_t n, const A &alpha,
	                 const T *x, size_t incx,
	                 const T *y, size_t incy, T *a, size_t lda);

### Arguments:

=m=
	Number of rows of `A`.
=n=
	Number of columns of `A`.
=alpha=
	The scale factor that is applied to the rank 1 product.
=x=
	Pointer to the first element of column vector `x`.
=incx=
	Increment (stride) of the column vector `x`.
=y=
	Pointer to the first element of row vector `x`.
=incy=
	Increment (stride) of the row vector `y`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).



---
## ConjugateRank1Update<a name="TBLAS_ConjugateRank1Update" />

Performs a rank 1 update of a matrix (`A <-- alpha * x * conj(y') + A`.
No index or dimension checking is performed.

This function corresponds to BLAS functions `zgerc` and `cgerc`.

### Prototype

	template <class A, class T>
	void ConjugateRank1Update(size_t m, size_t n, const A &alpha,
	                          const T *x, size_t incx,
	                          const T *y, size_t incy, T *a, size_t lda);

### Arguments:

=m=
	Number of rows of `A`.
=n=
	Number of columns of `A`.
=alpha=
	The scale factor that is applied to the rank 1 product.
=x=
	Pointer to the first element of column vector `x`.
=incx=
	Increment (stride) of the column vector `x`.
=y=
	Pointer to the first element of conjugated row vector `x`.
=incy=
	Increment (stride) of the conjugated row vector `y`.
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).



---
## MultMM<a name="TBLAS_MultMM" />

Performs a matrix-matrix multiply (`y <-- alpha * op(A) * op(B) + beta * C`, where `op` is specified by template parameters `transa` and `transb`).
No index or dimension checking is performed.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameters correspond to the classical `BLAS` parameters `transa` and `transb` and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::MultMM<'N','N'>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

This function corresponds to BLAS functions `zgemm`, `dgemm`, `cgemm`, and `sgemm`.

### Prototype

	template <char transa, char transb>
	struct MultMM{
		template <class A, class B, class T>
		MultMM(size_t m, size_t n, size_t k, const A &alpha, const T *a, size_t lda,
		       const T *b, size_t ldb, const B &beta, T *c, size_t ldc);
	};

### Arguments:

=transa=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `transa=N`, `op(A) = A`.
	* If `transa=T`, `op(A) = A'` (transpose of A)
	* If `transa=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=transb=
	Specifies the transpose operation used for `B`.
	This is a template parameter which must be explicitly specified.
	* If `transb=N`, `op(B) = B`.
	* If `transb=T`, `op(B) = B'` (transpose of B)
	* If `transb=C`, `op(B) = conj(B')` (conjugate transpose of B)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=m=
	Number of rows of `op(A)` or the result `C` (i.e. if `trans=T`, `m` is the number of columns of `A`).
=n=
	Number of columns of `op(B)` or the result `C`.
=alpha=
	The scale factor applied to `op(A) * op(B)` (for a simple matrix-matrix multiply, set this to 1).
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=b=
	Pointer to the first element of matrix `B`.
=ldb=
	Leading dimension of matrix `B`.
=beta=
	The scale factor applied to `C` before adding the product (for a simple matrix-matrix multiply, set this to 0).
=C=
	Pointer to the first element of destination matrix `C`.
=ldc=
	Leading dimension of destination matrix `C`.



---
## SolveTrM<a name="TBLAS_SolveTrM" />

Performs a matrix-matrix solve (solve for `X` in `op(A) * X = alpha * B` or `X * op(A) = alpha * B`, where `op` is specified by template parameter `transa`).
No index or dimension checking is performed.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameters correspond to the classical `BLAS` parameters `side`, `uplo`, `transa`, and `diag`, and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::SolveTrM<'L','U','N','N'>(m, n, alpha, a, lda, b, ldb);

This function corresponds to BLAS functions `ztrsm`, `dtrsm`, `ctrsm`, and `strsm`.

### Prototype

	template <char side, char uplo, char transa, char diag>
	struct SolveTrM{
		template <class TA, class T>
		SolveTrM(size_t m, size_t n, const TA &alpha, const T *a, size_t lda, T *b, size_t ldb);
	};

### Arguments:

=side=
	Whether `op(A)` appears on the left or right of `X`.
	This is a template parameter which must be explicitly specified.
	* If `side=L`, solves `op(A) * X = alpha * B`.
	* If `side=R`, solves `X * op(A) = alpha * B`.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=uplo=
	Specifies whether `A` is upper or lower triangular.
	This is a template parameter which must be explicitly specified.
	* If `uplo=U`, `A` is an upper triangular matrix.
	* If `uplo=L`, `A` is a lower triangular matrix.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=transa=
	Specifies the transpose operation used for `A`.
	This is a template parameter which must be explicitly specified.
	* If `transa=N`, `op(A) = A`.
	* If `transa=T`, `op(A) = A'` (transpose of A)
	* If `transa=C`, `op(A) = conj(A')` (conjugate transpose of A)
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these three uppercase characters.
=diag=
	Whether the diagonal of `A` is all ones (unit triangular).
	This is a template parameter which must be explicitly specified.
	* If `diag=U`, `A` is assumed to be unit triangular.
	* If `diag=N`, `A` is not assumed to be unit triangular.
	
	Exact checking may is not performed for this parameter, but it should always be specified as one of these two uppercase characters.
=m=
	Number of rows of `B`.
=n=
	Number of columns of `B`.
=alpha=
	The scale factor applied to result (for a simple triangular matrix solve, set this to 1).
=a=
	Pointer to the first element of matrix `A`.
=lda=
	Leading dimension of matrix `A` (typically the number of rows of `A` unless `A` is a submatrix of a larger matrix).
=b=
	Pointer to the first element of result matrix `B`.
=ldb=
	Leading dimension of result matrix `B`.




---
## Conjugate<a name="TBLAS_Conjugate" />

Conjugates a vector (`x <-- conj(x)`).

This function corresponds to LAPACK functions `zlacgv` and `clacgv`.

### Prototype

	template <class T>
	void Conjugate(size_t n, T *x, size_t incx);

### Arguments:

=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=incx=
	Increment (stride) of the vector.


---
## CopyMatrix<a name="TBLAS_CopyMatrix" />

Copies all or part of one matrix to another (`dst <-- src`).

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameter corresponds to the original `LAPACK` parameter `uplo` and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::CopyMatrix<'F'>(m, n, src, ldsrc, dst, lddst);

This function corresponds to LAPACK functions `zlacpy`, `dlacpy`, `clacpy`, and `slacpy`.

### Prototype

	template <char uplo>
	struct CopyMatrix{
		template <class T>
		CopyMatrix(size_t m, size_t n, const T *src, size_t ldsrc, T *dst, size_t lddst);
	};

### Arguments:

=uplo=
	Specifies which part of the matrices to copy.
	This is a template parameter which must be explicitly specified.
	* If `uplo=U`, only the upper triangle/trapezoid (including diagonal) is copied.
	* If `uplo=L`, only the lower triangle/trapezoid (including diagonal) is copied.
	* Any other value will cause the full matrix to be copied.
=m=
	Number of rows of the matrices.
=n=
	Number of columns of the matrices.
=src=
	Pointer to the first element of the source matrix.
=ldsrc=
	Leading dimension of the source matrix (typically the number of rows of `src` unless `src` is a submatrix of a larger matrix).
=dst=
	Pointer to the first element of the destination matrix.
=lddst=
	Leading dimension of the destination matrix.


---
## SetMatrix<a name="TBLAS_SetMatrix" />

Initializes the diagonal and off diagonal elements of all or part of a matrix.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameter corresponds to the original `LAPACK` parameter `uplo` and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::SetMatrix<'L'>(m, n, 0., 1., a, lda);

This function corresponds to LAPACK functions `zlaset`, `dlaset`, `claset`, and `slaset`.

### Prototype

	template <char uplo>
	struct SetMatrix{
		template <class TA, class TB, class T>
		SetMatrix(size_t m, size_t n, const TA &offdiag, const TB &diag, T *a, size_t lda);
	}

### Arguments:

=uplo=
	Specifies which part of the matrices to set.
	This is a template parameter which must be explicitly specified.
	* If `uplo=U`, only the upper triangle/trapezoid (including diagonal) is set.
	* If `uplo=L`, only the lower triangle/trapezoid (including diagonal) is set.
	* Any other value will cause the full matrix to be copied.
=m=
	Number of rows of the matrix.
=n=
	Number of columns of the matrix.
=a=
	Pointer to the first element of the matrix.
=lda=
	Leading dimension of the matrix (typically the number of rows of `a` unless `a` is a submatrix of a larger matrix).



---
## RandomVector<a name="TBLAS_RandomVector" />

Fills a vector with random complex numbers.

This function is doubly templated and is, in fact, implemented as the constructor of a `struct`.
The outer template parameter corresponds to the original `LAPACK` parameter `idist` and must be specified.
Here is an example of how the function should be called:

	RNP::TBLAS::RandomVector<2>(n, x);

This function corresponds to LAPACK functions `zlarnv` and `clarnv`.

### Prototype

	template <int dist>
	struct RandomVector{
		template <class T>
		RandomVector(size_t n, T *x, int iseed[4] = NULL);
	};

### Arguments:

=dist=
	Specifies the distribution of the random numbers.
	This is a template parameter which must be explicitly specified.
	* If `dist=1`, real and imaginary parts each uniformly selected from (0,1).
	* If `dist=2`, real and imaginary parts each uniformly selected from (-1,1).
	* If `dist=3`, real and imaginary parts each from standard normal distribution.
	* If `dist=4`, uniformly distributed on the unit disc `abs(z) < 1`.
	* If `dist=5`, uniformly distributed on the unit circle `abs(z) = 1`.
=n=
	Length of the vector.
=x=
	Pointer to the first element of the vector.
=iseed=
	The seed of the random number generator.
	Each of the 4 elements must be between 0 and 4095 (inclusive), and `iseed[3]` must be odd.
	The seed is updated on exit. If no seed is provided (by passing in `NULL`), an internal seed is used and updated.

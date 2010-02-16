% RNP Documentation -- Sparse
% Victor Liu (vkl@stanford.edu)
% January 16, 2010
<style type="text/css">
@import url(rnp.css);
</style>

# RNP::Sparse

The Sparse module contains sparse matrix structures, multiplication routines, and optional IO routines.
All functions and data structures are in namespace `RNP::Sparse`.

If the preprocessor flag `RNP_SPARSE_USE_IO`, then the IO functions are defined, and `IO.h` is required.

## Contents

### Sparse matrix structures
* [TCCSMatrix](#Sparse_TCCSMatrix)
* [TCRSMatrix](#Sparse_TCRSMatrix)

### Sparse matrix utility routines
* [MultMV](#Sparse_MultMV)
* [PrintSparseMatrix](#Sparse_PrintSparseMatrix)

---
## TCCSMatrix<a name="Sparse_TCCSMatrix" />

A compressed column storage sparse matrix, templated on the value type it contains.
This structure is not meant to be used directly; it should be constructed once from an entry-map, and the only subsequent operations it may perform are multiplying vectors.

### Constructor

There is only one constructor of any significance:

	template <class T>
	typedef std::pair<size_t,size_t> index_t;
	typedef std::map<index_t, T, index_comp_t> entry_map_t;
	TCCSMatrix(size_t N, const entry_map_t &entries);

This constructs an NxN matrix from a set of index-value pairs.
Each index is a pair of 0-based (row,column) indexes.
It is assumed that the matrix contains at least one nonzero element per column (or else the constructor fails silently).


---
## TCRSMatrix<a name="Sparse_TCRSMatrix" />

A compressed row storage sparse matrix, templated on the value type it contains.
This structure is not meant to be used directly; it should be constructed once from an entry-map, and the only subsequent operations it may perform are multiplying vectors.

### Constructor

There is only one constructor of any significance:

	template <class T>
	typedef std::pair<size_t,size_t> index_t;
	typedef std::map<index_t, T, index_comp_t> entry_map_t;
	TCRSMatrix(size_t N, const entry_map_t &entries);

This constructs an NxN matrix from a set of index-value pairs.
Each index is a pair of 0-based (row,column) indexes.
It is assumed that the matrix contains at least one nonzero element per row (or else the constructor fails silently).



---
## MultMV<a name="Sparse_MultMV" />

Performs a matrix-vector multiply (actually, an Axpy-like operation, `y <-- alpha * A * x + beta * y).
There is one overload defined for each type of matrix structure defined above.

### Prototype

	template <class T>
	void MultMV(
		const MATRIX_TYPE<T> &A, const T *X, T *Y,
		const T &scale_AX = T(1),
		const T &scale_Y = T(0)
		);
	
	The type `MATRIX_TYPE` is one of `TCCSMatrix` or `TCRSMatrix`.

### Arguments:

=A=
	The matrix `A`.
=X=
	Pointer to the first element of vector `x`.
	It is assumed that the elements are stored contiguously (`incx=1`)
=Y=
	Pointer to the first element of the result vector `y`.
	It is assumed that the elements are stored contiguously (`incy=1`)
=scale_AX=
	The scale factor for the product (`alpha` in the above equation).
=scale_Y=
	The scale factor for the original vector before adding the product (`beta` in the above equation).



---
## PrintSparseMatrix<a name="Sparse_PrintSparseMatrix" />

Outputs a sparse matrix.
Returns the output stream passed in.

This family of functions is only defined if `RNP_SPARSE_USE_IO` is defined.

### Prototype

	template <class T>
	void PrintSparseMatrix(const MATRIX_TYPE<T> &A);
	
	The type `MATRIX_TYPE` is one of `TCCSMatrix` or `TCRSMatrix`.

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



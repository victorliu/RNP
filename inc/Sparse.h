#ifndef _RNP_SPARSE_H_
#define _RNP_SPARSE_H_

#include <memory>
#include <map>
#include <complex>

// Preprocessor flags:
//   RNP_SPARSE_USE_IO

namespace RNP{
namespace Sparse{

template <typename T>
struct TCCSMatrix{
	int m,n;
	int *colptr; // points to start of columns in values
	int *rowind;
	T *values;

	TCCSMatrix():m(0),n(0),colptr(NULL),rowind(NULL),values(NULL){}
	TCCSMatrix(size_t N, size_t nnz):m(N),n(N),colptr(NULL),rowind(NULL),values(NULL){
		values = new T[nnz];
		colptr = new int[n+1];
		rowind = new int[nnz];
	}
	TCCSMatrix(size_t M, size_t N, size_t nnz):m(M),n(N),colptr(NULL),rowind(NULL),values(NULL){
		values = new T[nnz];
		colptr = new int[n+1];
		rowind = new int[nnz];
	}
	TCCSMatrix(const TCCSMatrix &M):m(M.m),n(M.n),colptr(NULL),rowind(NULL),values(NULL){
		int nnz = M.colptr[M.n];
		values = new T[nnz];
		colptr = new int[n+1];
		rowind = new int[nnz];
		
		std::uninitialized_copy(M.colptr, M.colptr+n+1, colptr);
		std::uninitialized_copy(M.rowind, M.rowind+nnz, rowind);
		std::uninitialized_copy(M.values, M.values+nnz, values);
	}
	
	typedef std::pair<size_t,size_t> index_t;
	struct index_comp_t{ // sort by columns first, then rows
		bool operator()(const index_t &a, const index_t &b) const{
			if(a.second < b.second){ return true; }
			else if(a.second > b.second){ return false; }
			else{ return a.first < b.first; }
		}
	};
	typedef std::map<index_t,T,index_comp_t> entry_map_t;
	// Assumes that entries contains at least one element per column.
	// Assumes all indexes in entries is consistent with the r and c given.
	// Assumes that flags is set consistent with entries.
	TCCSMatrix(size_t M, size_t N, const entry_map_t &entries):m(M),n(N){
		size_t nnz = entries.size();
		values = new T[nnz];
		colptr = new int[n+1];
		rowind = new int[nnz];
		
		size_t ip = 0;
		int prevcol = 0;
		colptr[0] = 0;
		for(typename entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			int col = i->first.second;
			rowind[ip] = i->first.first;
			values[ip] = i->second;
			
			++ip;
			if(prevcol != col){
				prevcol = col;
				colptr[col] = ip-1;
			}
		}
		colptr[n] = nnz; // do this at the end in case entries was bad, at least this is correct
	}
	TCCSMatrix(size_t N, const entry_map_t &entries):m(N),n(N){
		size_t nnz = entries.size();
		values = new T[nnz];
		colptr = new int[n+1];
		rowind = new int[nnz];
		
		size_t ip = 0;
		int prevcol = 0;
		colptr[0] = 0;
		for(typename entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			int col = i->first.second;
			rowind[ip] = i->first.first;
			values[ip] = i->second;
			
			++ip;
			if(prevcol != col){
				prevcol = col;
				colptr[col] = ip-1;
			}
		}
		colptr[n] = nnz; // do this at the end in case entries was bad, at least this is correct
	}
	TCCSMatrix& operator=(const TCCSMatrix &M){
		if(NULL != rowind){ delete [] rowind; }
		if(NULL != values){ delete [] values; }
		if(NULL != colptr){ delete [] colptr; }

		int n = M.n;
		colptr = new int[n+1];
		std::uninitialized_copy(M.colptr, M.colptr+n+1, colptr);

		int nnz = M.colptr[n];
		values = new T[nnz];
		rowind = new int[nnz];
		
		std::uninitialized_copy(M.rowind, M.rowind+nnz, rowind);
		std::uninitialized_copy(M.values, M.values+nnz, values);

		return *this;
	}
	virtual ~TCCSMatrix(){
		if(NULL != values){ delete [] values; }
		if(NULL != rowind){ delete [] rowind; }
		if(NULL != colptr){ delete [] colptr; }
	}
};

template <typename T>
struct TCRSMatrix{
	int m,n;
	int *colind;
	int *rowptr;
	T *values;

	TCRSMatrix():m(0),n(0),colind(NULL),rowptr(NULL),values(NULL){}
	TCRSMatrix(size_t M, size_t N, size_t nnz):m(M),n(N),colind(NULL),rowptr(NULL),values(NULL){
		rowptr = new int[m+1];
		values = new T[nnz];
		colind = new int[nnz];
	}
	TCRSMatrix(size_t N, size_t nnz):m(N),n(N),colind(NULL),rowptr(NULL),values(NULL){
		rowptr = new int[m+1];
		values = new T[nnz];
		colind = new int[nnz];
	}
	TCRSMatrix(const TCRSMatrix &M):m(M.m),n(M.n),colind(NULL),rowptr(NULL),values(NULL){
		int nnz = M.rowptr[m];
		values = new T[nnz];
		colind = new int[nnz];
		rowptr = new int[m+1];
		
		std::uninitialized_copy(M.colind, M.colind+nnz, colind);
		std::uninitialized_copy(M.rowptr, M.rowptr+m+1, rowptr);
		std::uninitialized_copy(M.values, M.values+nnz, values);
	}
	
	typedef std::pair<size_t,size_t> index_t;
	struct index_comp_t{ // sort by rows first, then columns
		bool operator()(const index_t &a, const index_t &b) const{
			if(a.first < b.first){ return true; }
			else if(a.first > b.first){ return false; }
			else{ return a.second < b.second; }
		}
	};
	typedef std::map<index_t,T,index_comp_t> entry_map_t;
	// Assumes that entries contains at least one element per column.
	// Assumes all indexes in entries is consistent with the r and c given.
	// Assumes that flags is set consistent with entries.
	TCRSMatrix(size_t N, const entry_map_t &entries):m(N),n(N){
		size_t nnz = entries.size();
		values = new T[nnz];
		rowptr = new int[m+1];
		colind = new int[nnz];
		
		size_t ip = 0;
		int prevrow = 0;
		rowptr[0] = 0;
		for(typename entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			int row = i->first.first;
			colind[ip] = i->first.second;
			values[ip] = i->second;
			
			++ip;
			if(prevrow != row){
				prevrow = row;
				rowptr[row] = ip-1;
			}
		}
		rowptr[m] = nnz; // do this at the end in case entries was bad, at least this is correct
	}
	TCRSMatrix(size_t M, size_t N, const entry_map_t &entries):m(M),n(N){
		size_t nnz = entries.size();
		values = new T[nnz];
		rowptr = new int[m+1];
		colind = new int[nnz];
		
		size_t ip = 0;
		int prevrow = 0;
		rowptr[0] = 0;
		for(typename entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			int row = i->first.first;
			colind[ip] = i->first.second;
			values[ip] = i->second;
			
			++ip;
			if(prevrow != row){
				prevrow = row;
				rowptr[row] = ip-1;
			}
		}
		rowptr[m] = nnz; // do this at the end in case entries was bad, at least this is correct
	}
	TCRSMatrix& operator=(const TCRSMatrix &M){
		if(NULL != colind){ delete [] colind; }
		if(NULL != values){ delete [] values; }
		if(NULL != rowptr){ delete [] rowptr; }

		int n = M.n;
		rowptr = new int[n+1];
		std::uninitialized_copy(M.rowptr, M.rowptr+n+1, rowptr);

		int nnz = M.rowptr[n];
		values = new T[nnz];
		colind = new int[nnz];
		
		std::uninitialized_copy(M.colind, M.colind+nnz, colind);
		std::uninitialized_copy(M.values, M.values+nnz, values);

		return *this;
	}
	virtual ~TCRSMatrix(){
		if(NULL != values){ delete [] values; }
		if(NULL != colind){ delete [] colind; }
		if(NULL != rowptr){ delete [] rowptr; }
	}
};

// Classes to differentiate complex types
template <class T>
struct _RealOrComplexChooser{
	typedef T value_type;
	typedef T real_type;

	inline static value_type _conj(const value_type &v){ return v; }
};
template <class T>
struct _RealOrComplexChooser<std::complex<T> >{
	typedef std::complex<T> value_type;
	typedef typename std::complex<T>::value_type real_type;

	inline static value_type _conj(const value_type &v){ using namespace std; return conj(v); }
};

template <char trans='N'>
struct MultMV{
	template <class T>
	MultMV(
		const TCCSMatrix<T> &A, const T *X, T *Y,
		const T &scale_AX = T(1),
		const T &scale_Y = T(0)
	){
		if('T' == trans){
			for(int i = 0; i < A.n; i++){
				T sum = 0;
				for(int jp = A.colptr[i]; jp < A.colptr[i+1]; jp++){
					int j = A.rowind[jp];
					sum += scale_AX*X[j]*A.values[jp];
				}
				Y[i] = scale_Y*Y[i] + scale_AX*sum;
			}
		}else if('C' == trans){
			for(int i = 0; i < A.n; i++){
				T sum = 0;
				for(int jp = A.colptr[i]; jp < A.colptr[i+1]; jp++){
					int j = A.rowind[jp];
					sum += scale_AX*X[j]*_RealOrComplexChooser<T>::_conj(A.values[jp]);
				}
				Y[i] = scale_Y*Y[i] + scale_AX*sum;
			}
		}else{
			for(int j = 0; j < A.n; j++){
				Y[j] *= scale_Y;
			}
			for(int j = 0; j < A.n; j++){
				for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
					int i = A.rowind[ip];
					Y[i] += scale_AX*X[j]*A.values[ip];
				}
			}
		}
	}
	template <class T>
	MultMV(
		const TCRSMatrix<T> &A, const T *X, T *Y,
		const T &scale_AX = T(1),
		const T &scale_Y = T(0)
	){
		if('T' == trans){
			for(int j = 0; j < A.n; j++){
				Y[j] *= scale_Y;
			}
			for(int j = 0; j < A.m; j++){
				for(int ip = A.rowptr[j]; ip < A.rowptr[j+1]; ip++){
					int i = A.colind[ip];
					Y[i] += scale_AX*X[j]*A.values[ip];
				}
			}
		}else if('C' == trans){
			for(int j = 0; j < A.n; j++){
				Y[j] *= scale_Y;
			}
			for(int j = 0; j < A.m; j++){
				for(int ip = A.rowptr[j]; ip < A.rowptr[j+1]; ip++){
					int i = A.colind[ip];
					Y[i] += scale_AX*X[j]*_RealOrComplexChooser<T>::_conj(A.values[ip]);
				}
			}
		}else{
			for(int i = 0; i < A.m; i++){
				T sum = 0;
				for(int jp = A.rowptr[i]; jp < A.rowptr[i+1]; jp++){
					int j = A.colind[jp];
					sum += scale_AX*X[j]*A.values[jp];
				}
				Y[i] = scale_Y*Y[i] + scale_AX*sum;
			}
		}
	}
};

#if defined(RNP_SPARSE_USE_IO)

#include "IO.h"

template <class T>
std::ostream& PrintSparseMatrix(const TCCSMatrix<T> &A, std::ostream &os = std::cout){
	#ifdef RNP_OUTPUT_MATHEMATICA
	os << "SparseArray[{";
	for(int j = 0; j < (int)A.n; j++){
		for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
			int i = A.rowind[ip];
			os << "{" << i+1 << ", " << j+1 << "} -> ";
			RNP::IO::Print(A.values[ip], os);
			os << ", ";
		}
	}
	os << "{_, _} -> 0}]";
	#endif
	#ifdef RNP_OUTPUT_MATLAB
	os << "spconvert([";
	for(int j = 0; j < (int)A.n; j++){
		for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
			int i = A.rowind[ip];
			os << i+1 << "\t" << j+1 << "\t";
			RNP::IO::Print(A.values[ip], os);
			os << std::endl;
		}
	}
	os << "])";
	#endif
	return os;
}

template <class T>
std::ostream& PrintSparseMatrix(const TCRSMatrix<T> &A, std::ostream &os = std::cout){
	#ifdef RNP_OUTPUT_MATHEMATICA
	os << "SparseArray[{";
	for(int i = 0; i < (int)A.m; i++){
		for(int jp = A.rowptr[i]; jp < A.rowptr[i+1]; jp++){
			int j = A.colind[jp];
			os << "{" << i+1 << ", " << j+1 << "} -> ";
			RNP::IO::Print(A.values[jp], os);
			os << ", ";
		}
	}
	os << "{_, _} -> 0}]";
	#endif
	#ifdef RNP_OUTPUT_MATLAB
	os << "spconvert([";
	for(int i = 0; i < (int)A.m; i++){
		for(int jp = A.colptr[i]; jp < A.colptr[i+1]; jp++){
			int j = A.rowind[jp];
			os << i+1 << "\t" << j+1 << "\t";
			RNP::IO::Print(A.values[jp], os);
			os << std::endl;
		}
	}
	os << "])";
	#endif
	return os;
}

#endif // defined(RNP_SPARSE_USE_IO)

}; // namespace Sparse
}; // namespace RNP

#endif // _RNP_SPARSE_H_

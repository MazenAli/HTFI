/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef LAWA_MATRIX_TOOLS_H
#define LAWA_MATRIX_TOOLS_H 1

#include <flens/flens.h>
#include <lawa/lawa.h>
#include <htucker/htucker.h>


using namespace flens;
using namespace htucker;
namespace lawa{



//template <typename T> 
//GeMatrix<FullStorage<T,ColMajor> >
//operator*(SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > & mat, GeMatrix<FullStorage<T,ColMajor> > & gemat){
//	assert(mat.numCols() == gemat.numRows());
//	GeMatrix<FullStorage<T,ColMajor> > ret(mat.numRows(),gemat.numCols());
//	flens::extensions::CRS<T,flens::CRS_General> eng = mat.engine();
//	T val = 0;
//	for(int i = 1; i <= mat.numRows(); ++i){
//		for(int j = 1; j<= gemat.numCols(); ++j){
//			val = 0;
//			for(int k = eng.rows(i); k < eng.rows(i+1); ++k){
//				val += eng.values(k)*gemat(eng.columns(k),j);
//			}
//			ret(i,j) = val;
//		}
//	}
//	return ret;
//};

template <typename T>
T getElement(SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > mat, int row, int col){
	DenseVector<Array<T> > v = (mat.engine()).values;
	DenseVector<Array<int> > r = (mat.engine()).rows;
	DenseVector<Array<int> > c = (mat.engine()).columns;
	int end = r.lastIndex();
	int searchend;
	if(row < end){
		searchend = r(row+1) - 1;
	} else {
		searchend = c.lastIndex();
	}

	for(int i = r(row); i <= searchend; ++i){
		if(c(i) == col){
			return v(i);
		}
	}
	return 0.0;
}

template <typename T>
SparseGeMatrix<flens::extensions::CRS<T> > 
kronecker_product(const SparseGeMatrix<flens::extensions::CRS<T> > &A,const SparseGeMatrix<flens::extensions::CRS<T> > &B)
{
	const flens::extensions::CRS<T> A_crs=A.engine(), B_crs=B.engine();
	int m=A.numRows(), n=A.numCols(),s=B.numRows(), t=B.numCols();
	SparseGeMatrix<flens::extensions::CRS<T> >C(m*s,n*t);
	
	for (int i=1; i<=m; ++i) {
		for (int j=A_crs.rows(i); j<=A_crs.rows(i+1)-1; ++j) {
			for (int k=1; k<=s; ++k) {
				for (int l=B_crs.rows(k); l<=B_crs.rows(k+1)-1; ++l) {
					C((i-1)*s+k, (A_crs.columns(j)-1)*t+B_crs.columns(l))=A_crs.values(j)*B_crs.values(l);
				}
			}
		}
	}
	C.finalize();
	return C;
}

//template <typename T>
//ostream& operator<<(ostream& Stream, SparseGeMatrix<flens::extensions::CRS<T> > &B)
//{
//	Stream << B.numRows() << " x " << B.numCols() << ": " <<  endl;
//	const DenseVector<Array<T> > & ref = (B.engine().values);
//	const DenseVector<Array<T> > & col = (B.engine().columns);
//	const DenseVector<Array<T> > & refrows = B.engine().rows ;
//	int j = 1;
//	for(int i = ref.firstIndex(); i <= ref.lastIndex(); ++i){
//		
//		if(j < refrows.length() && refrows(j+1) == i){
//			j = j + 1;
//		}
//		Stream << "(" << j << ", " << col(i) << ") = " << ref(i) << endl;
//	}
//	
//	return Stream;
//};

template <typename T>
GeMatrix<FullStorage<T,ColMajor> > 
kronecker_product(const GeMatrix<FullStorage<T,ColMajor> >  &A,const GeMatrix<FullStorage<T,ColMajor> >  &B)
{
	GeMatrix<FullStorage<T,ColMajor> >C(A.numRows()*B.numRows(),A.numCols()*B.numCols());
	
	for(int i1 = A.firstRow(); i1 <= A.lastRow(); i1++){
		for(int i2 = B.firstRow(); i2 <= B.lastRow(); i2++){
			for(int i3 = A.firstCol(); i3 <= A.lastCol(); i3++){
				for(int i4 = B.firstCol(); i4 <= B.lastCol(); i4++){
					C((i1 - A.firstRow())*B.numRows() + i2, (i3-A.firstCol())*B.numCols() + i4) = A(i1,i3)*B(i2,i4);
				}
			}
		}
	}
	return C;
};

template <typename T>
SparseGeMatrix<flens::extensions::CRS<T> > 
subMatrix(const SparseGeMatrix<flens::extensions::CRS<T> > &A, const Range<int> r1, const Range<int> r2)
{
	const flens::extensions::CRS<T> A_crs=A.engine();
	SparseGeMatrix<flens::extensions::CRS<T> > C(r1.length(),r2.length());
	int row_index = 1;
	for(int i = A_crs.values.firstIndex(); i<= A_crs.values.lastIndex(); ++i){
		while(!(A_crs.rows(row_index) <= i && A_crs.rows(row_index+1) > i)){
			row_index ++;
		}
		if(r1.firstIndex() <= row_index && r1.lastIndex() >= row_index && r2.firstIndex() <= A_crs.columns(i) && A_crs.columns(i) <= r2.lastIndex()){
			C(row_index-r1.firstIndex() + 1, A_crs.columns(i) - r2.firstIndex() + 1) = A_crs.values(i);
		}
	}
	C.finalize();
	return C;
}

template <typename T>
DenseVector<Array<T> >
getDiagonal(const SparseGeMatrix<flens::extensions::CRS<T> > &A){
	const flens::extensions::CRS<T> A_crs = A.engine();
	assert(A_crs.numRows() == A_crs.numCols());
	DenseVector<Array<T> > ret(A_crs.numRows());
	cout << A_crs.rows.firstIndex() << endl;
	
	for(int i = A_crs.rows.firstIndex(); i < A_crs.rows.lastIndex(); ++i){
		int search = A_crs.rows(i);
		while(search < A_crs.rows(i+1) && A_crs.columns(search) < i){
			search ++;
		}
		if(A_crs.columns(search) == i){
			ret(i) = A_crs.values(search);
		} else {
			ret(i) = 0.0;
		}
		//cout << "i = " << i << "  " <<  A_crs.rows(i) << endl;
	}
	return ret;
}

template <typename T>
SparseGeMatrix<flens::extensions::CRS<T> >
Unitmatrix(int size){
	SparseGeMatrix<flens::extensions::CRS<T> > ret(size,size);
	for(int i = 1; i <= size; ++i){
		ret(i,i) = 1;
	}
	ret.finalize();
	return ret;
}

template <typename T>
SparseGeMatrix<flens::extensions::CRS<T> >
operator+ (const SparseGeMatrix<flens::extensions::CRS<T> > &A,const SparseGeMatrix<flens::extensions::CRS<T> > &B){
	DenseVector<Array<T> > v = (A.engine()).values;
	DenseVector<Array<int> > r = (A.engine()).rows;
	DenseVector<Array<int> > c = (A.engine()).columns;
	SparseGeMatrix<flens::extensions::CRS<T> > C(A.numRows(),A.numCols());
	int rcount = r.firstIndex();
	for(int i = v.firstIndex(); i <= v.lastIndex(); ++i){
		if(r(rcount + 1) <= i){
			rcount++;
		}
		C(rcount,c(i)) = v(i);
	}
	v = (B.engine()).values;
	r = (B.engine()).rows;
	c = (B.engine()).columns;
	rcount = r.firstIndex();
	for(int i = v.firstIndex(); i <= v.lastIndex(); ++i){
		if(r(rcount + 1) <= i){
			rcount++;
		}
		C(rcount,c(i)) += v(i);
	}
	C.finalize();
	return C;
}

template <typename T>
SparseGeMatrix<flens::extensions::CRS<T> >
operator* (const double constval ,const SparseGeMatrix<flens::extensions::CRS<T> > &B){
	SparseGeMatrix<flens::extensions::CRS<T> > C(B.numRows(),B.numCols());
	C = B;
	for(int i = C.engine().values.firstIndex(); i <= C.engine().values.lastIndex(); ++i){
		C.engine().values(i) *= constval;
	}
	
	return C;
}

template <typename T>
T LCDet(GeMatrix<FullStorage<T,ColMajor> > mat){
	assert(mat.numCols() == mat.numRows());
	GeMatrix<FullStorage<T,ColMajor> > last,lastlast,akt;
	akt = mat;
	for(int i = mat.numRows()-1; i>=1; --i){
		lastlast = last;
		//cout << "lastlast = " << lastlast << endl;
		last = akt;
		akt = GeMatrix<FullStorage<T,ColMajor> >(i,i);
		for(int j = 1; j <= i ; ++j){
			for(int k = 1; k <= i ; ++k){
				akt(j,k) = last(j,k)*last(j+1,k+1) - last(j+1,k)*last(j,k+1);
				if(i < mat.numRows() -1){
						akt(j,k) /= lastlast(j+1,k+1);
				}
			}
		}
	}
	return akt(1,1);
}

template <typename T>
T Det(typename TrMatrix<FullStorage<T,ColMajor> >::View mat){
	assert(mat.lastCol() - mat.firstCol() == mat.lastRow()- mat.firstRow());
	T prod = 1.0;
	for(int i = mat.firstCol(); i<= mat.lastCol(); ++i){
		prod *= mat(i,i);
	}
	return prod;
}

template <typename T>
T Det(GeMatrix<FullStorage<T,ColMajor> > mat){
	assert(mat.numRows() == mat.numCols());
	if(mat.numRows() == 1){
		return mat(1,1);
	} else {
		GeMatrix<FullStorage<T,ColMajor> > tmp = mat;
		DenseVector<Array<int> > piv(mat.numRows());
		int count = 0;
		trf(tmp, piv);
		for(int i = piv.firstIndex(); i<= piv.lastIndex(); ++i){
			if(i != piv(i)) {
				count++;
			}
		}
		GeMatrix<FullStorage<double,ColMajor> >::TriangularView U = upper(tmp);
		
		count = (count %2)*(-2) + 1;

		return Det<double>(U)*(count+0.0);
	}
};


template <typename T>
T LCMinor(GeMatrix<FullStorage<T,ColMajor> > mat,int _row, int _col){
	assert(mat.numCols() == mat.numRows());
	GeMatrix<FullStorage<T,ColMajor> > akt;
	if(mat.numRows()>1){
		akt = GeMatrix<FullStorage<T,ColMajor> >(mat.numRows()-1,mat.numCols()-1);
		for(int i = 1; i <= mat.numRows()-1; ++i){
			for(int j = 1; j <= mat.numCols()-1; ++j){
				if(i < _row){
					if(j < _col){
						akt(i,j) = mat(i,j);
					}  else {
						akt(i,j) = mat(i,j+1);
					}
				} else {
					if(j < _col){
						akt(i,j) = mat(i+1,j);
					}  else {
						akt(i,j) = mat(i+1,j+1);
					}
				}
			
			}
		}
		return LCDet<T>(akt);
	} else {
		return T(1);
	}
}

//template <typename T>
//GeMatrix<FullStorage<T,ColMajor> > Inverse(GeMatrix<FullStorage<T,ColMajor> > mat){
//	T det = LCDet<T>(mat);
//	GeMatrix<FullStorage<T,ColMajor> > ret = mat;
//	for(int i = 1; i <= mat.numRows(); ++i){
//		for(int j = 1; j <= mat.numCols(); ++j){
//			ret(i,j) = LCMinor<T>(mat,j,i)/det;
//			if((i+j) % 2 != 0){
//				ret(i,j) *= -1;
//			}
//		}
//	}
//	return ret;
//}




template <typename T>
GeMatrix<FullStorage<T,ColMajor> > Inverse(GeMatrix<FullStorage<T,ColMajor> > mat){
	assert(mat.numRows() == mat.numCols());
	GeMatrix<FullStorage<T,ColMajor> > tmp = mat;
	DenseVector<Array<int> > piv(mat.numRows());
	trf(tmp,piv);
    tri(tmp,piv);
	return tmp;	
}

template <typename T>
GeMatrix<FullStorage<T,ColMajor> > Vec2Mat(DenseVector<Array<T> > vec,int m, int n){
	assert(vec.length() == m*n);
	GeMatrix<FullStorage<T,ColMajor> > ret(m,n);
	for(int i = 1;i <= m; i++){
		for(int j = 1; j <= n; j++){
			ret(i,j) = vec((j-1)*m+i);
		}
	}
	return ret;
};

template <typename T>
GeMatrix<FullStorage<T,ColMajor> > Vec2Mat(typename DenseVector<Array<T> >::View vec,int m, int n){
	assert(vec.length() == m*n);
	GeMatrix<FullStorage<T,ColMajor> > ret(m,n);
	for(int i = 1;i <= m; i++){
		for(int j = 1; j <= n; j++){
			ret(i,j) = vec((j-1)*m+i);
		}
	}
	return ret;
};

template <typename T>
DenseVector<Array<T> > Mat2Vec(GeMatrix<FullStorage<T,ColMajor> > mat){
	int rows = (mat.lastRow() - mat.firstRow() + 1);
	int cols = (mat.lastCol() - mat.firstCol() + 1);
	DenseVector<Array<T> > ret(rows*cols);
	for(int i = mat.firstRow(); i <= mat.lastRow(); ++i){
		for(int j = mat.firstCol(); j<= mat.lastCol(); ++j){

			ret((j-mat.firstCol())*(rows) + i) = mat(i,j);
		}
	}
 return ret;
};



//nicht schön dafür extrem selten:
template <typename T>
void
mm(Transpose transA,
   Transpose transB,
   T alpha,
   const GeMatrix<FullStorage<T, ColMajor> > &A,
   typename TrMatrix<FullStorage<T, ColMajor> >::ConstView &B,
   T beta,
   GeMatrix<FullStorage<T, ColMajor> > &C)
{

	if(B.upLo() == Upper){
	   GeMatrix<FullStorage<T, ColMajor> > D(B.lastRow() - B.firstRow() +1 ,B.lastCol() - B.firstCol() +1);
	   for(int i = B.firstRow(); i <= B.lastRow(); ++i){
		   for(int j = B.firstCol(); j <= B.lastCol(); ++j){
			   if(j >= i){
					D(i,j) = B(i,j);
			   } else {
					D(i,j) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,A,D,beta,C);
	} else {
		GeMatrix<FullStorage<T, ColMajor> > D(B.lastRow() - B.firstRow() +1 ,B.lastCol() - B.firstCol() +1);
	    for(int i = B.firstRow(); i <= B.lastRow(); ++i){
		   for(int j = B.firstCol(); j <= B.lastCol(); ++j){
			   //cout << i << " " << j << endl;	
			   if(j <= i){
				D(i,j) = B(i,j);
			   } else {
				D(i,j ) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,A,D,beta,C);
	}
};

template <typename T>
void
mm(Transpose transA,
   Transpose transB,
   T alpha,
   typename TrMatrix<flens::FullStorage<T, ColMajor> >::ConstView &A,
   const GeMatrix<FullStorage<T, ColMajor> > &B,
   T beta,
   GeMatrix<FullStorage<T, ColMajor> > &C)
{
   if(A.upLo() == Upper){
	   GeMatrix<FullStorage<T, ColMajor> > D(A.lastRow() - A.firstRow() +1 ,A.lastCol() - A.firstCol() +1);
	   for(int i = A.firstRow(); i <= A.lastRow(); ++i){
		   for(int j = A.firstCol(); j <= A.lastCol(); ++j){
			   if(j >= i){
					D(i,j) = A(i,j);
			   } else {
					D(i,j) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,D,B,beta,C);
	} else {
		GeMatrix<FullStorage<T, ColMajor> > D(A.lastRow() - A.firstRow() +1 ,A.lastCol() - A.firstCol() +1);
	    for(int i = A.firstRow(); i <= A.lastRow(); ++i){
		   for(int j = A.firstCol(); j <= A.lastCol(); ++j){
			   if(j <= i){
				D(i,j) = A(i,j);
			   } else {
				D(i,j ) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,D,B,beta,C);
	}
};


template <typename T>
GeMatrix<FullStorage<T,ColMajor> >
sparse2full(SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > &A){
	flens::extensions::CRS<T,flens::CRS_General> & cs = A.engine();
	GeMatrix<FullStorage<T,ColMajor> > B(A.numRows(),A.numCols());
	int j;
	for(int i = 1; i <= A.numRows(); ++i){
		int end;
		if(i == A.numRows()){
			end = cs.columns.length() + 1;
		} else {
			end = cs.rows(i+1);
		};
		j = cs.rows(i);
		while(j < end){
			B(i,cs.columns(j)) = cs.values(j);
			++j;
		};
		
	};

	return B;
};

template <typename T> 
GeMatrix<FullStorage<T,ColMajor> >
tovector(GeMatrix<FullStorage<T,ColMajor> > & A){
	GeMatrix<FullStorage<T,ColMajor> > B(A.numRows()*A.numCols(),1);
	for(int i = 1; i <= A.numCols(); ++i){
		B(_((i-1)*A.numRows()+1,i*A.numRows()),1) = A(_,i);
	}
	return B;
};


template <typename T>
void
mm(Transpose transA,
   Transpose transB,
   T alpha,
   typename TrMatrix<FullStorage<T, ColMajor> >::ConstView &A,
   typename TrMatrix<FullStorage<T, ColMajor> >::ConstView &B,
   T beta,
   GeMatrix<FullStorage<T, ColMajor> > &C)
{
   if(A.upLo() == Upper){
	   GeMatrix<FullStorage<T, ColMajor> > D(A.lastRow() - A.firstRow() +1 ,A.lastCol() - A.firstCol() +1);
	   for(int i = A.firstRow(); i <= A.lastRow(); ++i){
		   for(int j = A.firstCol(); j <= A.lastCol(); ++j){
			   if(j >= i){
					D(i,j) = A(i,j);
			   } else {
					D(i,j) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,D,B,beta,C);
	} else {
		GeMatrix<FullStorage<T, ColMajor> > D(A.lastRow() - A.firstRow() +1 ,A.lastCol() - A.firstCol() +1);
	    for(int i = A.firstRow(); i <= A.lastRow(); ++i){
		   for(int j = A.firstCol(); j <= A.lastCol(); ++j){
			   if(j <= i){
				D(i,j) = A(i,j);
			   } else {
				D(i,j ) = 0;
			   }
		   }
	   }
		mm(transA,transB,alpha,D,B,beta,C);
	}
};


} //namespace lawa

#endif // LAWA_MATRIX_TOOLS_H

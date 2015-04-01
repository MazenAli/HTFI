namespace htucker{

template <typename T, typename TensorFunction>
ROA<T, TensorFunction>::ROA(const double _eps):d(0),Atensor(1,1),U(1,1),V(1,1),SVD_epsilon(_eps){};

template <typename T, typename TensorFunction>
ROA<T, TensorFunction>::ROA(const tensor<TensorFunction> &_A, const DimensionIndex &_t, const DimensionIndex &_tbar, const int _dim, const double _eps):A(_A),t(_t),tbar(_tbar),d(_dim),Atensor(1,1),U(1,1),V(1,1),tcomplement(_t.getComplement(_dim)),SVD_epsilon(_eps){};


template <typename T, typename TensorFunction>
ROA<T, TensorFunction>::ROA(const ROA<T,TensorFunction> &copy):t(copy.getT()),tbar(copy.getTbar()),tcomplement(copy.getT().getComplement(copy.dim())),d(copy.dim()),A(copy.A),Atensor(copy.Atensor),U(copy.U),V(copy.V),SVD_epsilon(copy.SVD_epsilon){}



template <typename T, typename TensorFunction>
void 
ROA<T, TensorFunction>::addpivot(const DimensionIndex &pivot){
	DimensionIndex ausw;
	this->pivots.add(pivot);

	//std::cout << "roa.tcc ROA<T, SVD>::addpivot(DimensionIndex pivot): " << pivots << std::endl;
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(pivots.length(),pivots.length());
	if(pivots.length() > 1){
		mat(_(1,pivots.length()-1),_(1,pivots.length()-1)) = Atensor;
	}
	
	Atensor = mat;
	assert(pivot.length() == d);
	ausw = pivots[pivots.length() -1];
	
	for(int i = 0; i < pivots.length(); ++i){
		ausw.setValue(pivots[i],t);
		Atensor(i+1,pivots.length()) = A(ausw);
	}

	ausw =pivots[pivots.length() -1];

	for(int j = 0; j < pivots.length()-1; ++j){
		ausw.setValue(pivots[j],tcomplement);
		Atensor(pivots.length(),j+1) = A(ausw);
	}
	//std::cout << "roa.tcc addpivot: vor SVD: Atensor = " << Atensor << std::endl;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Atensave = Atensor;
	flens::svd(Atensave,s,U,V);
	inverse_s = flens::DenseVector<flens::Array<T> >(s.length());
	for(int i = s.firstIndex(); i<= s.lastIndex(); ++i){
		if(abs(s(i)) > SVD_epsilon){
			inverse_s(i) = 1/s(i);
		} 
	}

	//std::cout << "roa.tcc ROA<T, SVD>::addpivot(DimensionIndex pivot), svd: s = " << s << " U = " << U << " VT = " << V << std::endl;
};

template <typename T, typename TensorFunction >
T 
ROA<T, TensorFunction >::operator()(const DimensionIndex & vals) const{
	double eps = SVD_epsilon;
	if(pivots.length() == 0) return 0;
	assert(vals.length() == d);
	T sum = 0;
	
	flens::DenseVector<flens::Array<T> > sol(pivots.length()); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
	DimensionIndex idx = vals;
	
	for(int i = 0; i < pivots.length(); ++i){
		idx.setValue(pivots[i],t);
		sol(i+1) = A(idx);
	}
	
		//std::cout << "roa.tcc operator() Atensor = " << Atensor << std::endl;
	flens::DenseVector<flens::Array<T> > savesol(pivots.length());
	
	savesol = flens::transpose(U)*sol;
	//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,T(1),U,sol,T(0),savesol);
	sol = flens::DenseVector<flens::Array<T> >(savesol.length());
	for(int i = savesol.firstIndex(); i <= savesol.lastIndex(); ++i){
		sol(i) = inverse_s(i)*savesol(i);
	}
	savesol = flens::transpose(V)*sol;
	sol = savesol;
	idx = vals;
	for(int i = 0; i < sol.length(); ++i){
		idx.setValue(pivots[i],tcomplement);
		//std::cout << idx << " A(idx) = " << A(idx) << " sol(i + 1,1) = " << sol(i+1,1) << "  " << tcomplement <<   std::endl;
		sum+= A(idx)*sol(i + 1);
	}
	return sum;
};

template <typename T, typename TensorFunction>
flens::DenseVector<flens::Array<T> >
ROA<T, TensorFunction>::evaluate(const DimensionIndex & initialIndex, const DimensionIndex &_t, const int dim, const int _min, const int _max) const{
	double eps = SVD_epsilon;
	
	if(pivots.length() == 0){
		flens::DenseVector<flens::Array<T> > ret(_max - _min + 1);
		for(int i = ret.firstIndex(); i <= ret.lastIndex(); ++i){
			ret(i) = 0;
		}
		return ret;
	}
	
	bool in_t = false;
	
	for(int i = 0; i < _t.length() ; ++i){
		if(t[i] == dim){
			in_t = true;
			break;
		}
	}
	
	if(in_t){
		//Die variable Dimension ist ein Zeilenindex
		flens::DenseVector<flens::Array<T> > Aright(pivots.length());
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Aleft(_max - _min + 1, pivots.length());
		DimensionIndex idx = initialIndex;
		//std::cout << "pivots.length() = " << pivots.length() << "pivots = " << pivots << std::endl;
		for(int i = 0; i < pivots.length(); ++i){
			idx.setValue(pivots[i],_t);
			//std::cout << "Aright: " << idx << std::endl;
			Aright(i+1) = A(idx);
		}
		
		for(int i = 0; i< pivots.length(); ++i){
			idx = pivots[i];
			idx.setValue(initialIndex,_t);
			for(int j = _min; j<= _max; ++j){
				idx.setValue(dim,j);
				//std::cout << "Aleft " << (j - _min + 1) << " ," <<  i + 1 << " = " << idx <<  std::endl;
				Aleft(j - _min + 1, i + 1) = A(idx);
			}
		}

		//std::cout << "Aleft = " << Aleft << std::endl;
		//std::cout << "Aright = " << Aright << std::endl;
		//std::cout << "s = " << s << std::endl;
		flens::DenseVector<flens::Array<T> > res1(U.numCols());
		res1 = flens::transpose(U) * Aright;
		for(int i = 1; i <= res1.length();++i){
			if(abs(s(i)) > eps){
				res1(i) = res1(i)/s(i);
			} else {
				res1(i) = 0;
			}
		}

		flens::DenseVector<flens::Array<T> > res2(V.numRows());
		res2 = flens::transpose(V) * res1;
		res1 = Aleft * res2;

		return res1;
	} else {
		//Die variable Dimension ist ein Spaltenindex
		flens::DenseVector<flens::Array<T> > Aleft(pivots.length());
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Aright( pivots.length(),_max - _min + 1);
		DimensionIndex idx = initialIndex;
		for(int i = 0; i < pivots.length(); ++i){
			idx = pivots[i];
			idx.setValue(initialIndex,_t);
			Aleft(i+1) = A(idx);
		}

		for(int i = 0; i< pivots.length(); ++i){
			idx = initialIndex ;
			idx.setValue(pivots[i],_t);
			for(int j = _min; j<= _max; ++j){
				idx.setValue(dim,j);
				Aright(i + 1,j - _min + 1) = A(idx);
			}
		}

		flens::DenseVector<flens::Array<T> > res1(V.numCols());
		res1 = V * Aleft;
		
		for(int i = 1; i <= res1.length();++i){
			if(abs(s(i)) > eps){
				res1(i) = res1(i)/s(i);
			} else {
				res1(i) = 0;
			}
		}

		flens::DenseVector<flens::Array<T> > res2(U.numRows());
		res2 = U * res1;
		res1 = flens::DenseVector<flens::Array<T> >(Aright.numCols());
		res1 = flens::transpose(Aright) * res2;
	
		return res1;
	}
	
};

template <typename T, typename TensorFunction>
flens::DenseVector<flens::Array<T> >
ROA<T, TensorFunction>::evaluate(const DimensionIndex & initialIndex, const DimensionIndex & _t, const DimensionIndexList &list, const DimensionIndex & list_dims)const {
	// we assume that all Dimensions in list_dims (father_pivot_dimensions) are not in _t (i.e. they are column indices!)
	double eps = SVD_epsilon;
	if(pivots.length() == 0){
		flens::DenseVector<flens::Array<T> > ret(list.length());
		for(int i = ret.firstIndex(); i <= ret.lastIndex(); ++i){
			ret(i) = 0;
		}

		return ret;
	}

	flens::DenseVector<flens::Array<T> > Aleft(pivots.length());
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Aright(pivots.length(), list.length());
	DimensionIndex idx = initialIndex;
	for(int i = 0; i < pivots.length(); ++i){
		idx = pivots[i];
		idx.setValue(initialIndex,_t);
		Aleft(i + 1) = A(idx);
	}
	for(int i = 0; i <pivots.length(); ++i){
		idx = initialIndex;
		idx.setValue(pivots[i],_t);
		for(int j = 0; j< list.length(); ++j){
			idx.setValue(list[j],list_dims);
			Aright(i+1,j+1) = A(idx);
		}
	}

	flens::DenseVector<flens::Array<T> > res1(V.numCols());
	res1 = V * Aleft;

	for(int i = 1; i <= res1.length();++i){
			if(abs(s(i)) > eps){
				res1(i) = res1(i)/s(i);
			} else {
				res1(i) = 0;
			}
		}

	flens::DenseVector<flens::Array<T> > res2(U.numRows());
	res2 = U * res1;
	res1 = flens::DenseVector<flens::Array<T> >(Aright.numCols());
	res1 = flens::transpose(Aright) * res2;

	return res1;
};



template <typename T, typename TensorFunction>
DimensionIndex  
ROA<T, TensorFunction>::GreedyPivotSearch(const int lmax, const DimensionIndexList & fatherpivots, double &error) const{
	#ifdef _GREEDYDEBUG
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, d = " << this->d << std::endl;
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, t = " << t << std::endl;
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, tbar = " << tbar << std::endl;
    #endif
	  DimensionIndex tjtbar = t.join(tbar);
	#ifdef _GREEDYDEBUG 
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, tjtbar = " << tjtbar << std::endl;
	#endif
	  DimensionIndex tjtbarcomplement;
	  if(tjtbar.length() < d){
		 tjtbarcomplement = tjtbar.getComplement(d);
	  }
	  T compare;
    #ifdef _GREEDYDEBUG
	  std::cout << "fatherpivots = " << fatherpivots << std::endl;
	  std::cout << "tjtbarcomplement = " << tjtbarcomplement << std::endl;
    #endif
	  DimensionIndex index(d);
	  int count = 0;
	
	  index.setRandom(tjtbar,A.getmin(),A.getmax());
	  if(tjtbarcomplement.length() > 0){
		  index.setRandom(fatherpivots,tjtbarcomplement);
	  }
	  DimensionIndex saveIndex(d);
	#ifdef _GREEDYDEBUG
	  std::cout << "index = " << index << "A.minval = " << A.getmin() << " A.maxval = " << A.getmax() <<  std::endl;	 
	  std::cout << "A(index) = " << A(index) << std::endl;
	  std::cout << "this->operator()(index) = " << this->operator()(index) << std::endl;
	#endif
	  double startmaximizer = abs(A(index) - this->operator()(index));
	  double savemaximizer;

	  for(int i = 1; i<= 10000; i++){
		  saveIndex.setRandom(tjtbar,A.getmin(),A.getmax());
		  if(tjtbarcomplement.length() > 0){
			saveIndex.setRandom(fatherpivots,tjtbarcomplement);
		  }
		  savemaximizer = abs(A(saveIndex) - this->operator()(saveIndex));
		  if(savemaximizer > startmaximizer){
			index = saveIndex;
			startmaximizer = savemaximizer;
		  }
	  }
	#ifdef _GREEDYDEBUG
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, index = " << index << std::endl;
	#endif
	  compare = abs(A(index) - (*this)(index));
	#ifdef _GREEDYDEBUG	 
	  std::cout << "roa.tcc ROA<T>::GreedyPivotSearch, compare = " << compare << std::endl;
	#endif
	  flens::DenseVector<flens::Array<T> > values;
	  for(int l = 1; l <= lmax; ++l){
		 for(int i = 0; i < tjtbar.length(); ++i){
			 values = this->evaluate(index,t,tjtbar[i],A.getmin()[tjtbar[i]-1],A.getmax()[tjtbar[i]-1]);
	#ifdef _GREEDYDEBUG
			 std::cout << "index = " << index <<  std::endl;
				 //std::cout << "values: " << values << std::endl;
    #endif
			 if(A.vecEval()){
				 #ifdef _GREEDYDEBUG
					std::cout << "index = " << index << std::endl;
					std::cout << tjtbar[i] << std::endl;
				 #endif
				 flens::DenseVector<flens::Array<T> > tensoreval = A.vec(index,tjtbar[i]);
				 int pos;
				 #ifdef _GREEDYDEBUG
					std::cout << "here" << std::endl;
				 #endif
				 flens::dvmaxabs<T>(tensoreval - values,pos);
				 #ifdef _GREEDYDEBUG
					std::cout << "here" << std::endl;
					std::cout << "pos-1 + A.getmin()[tjtbar[i]-1] = " << pos-1 + A.getmin()[tjtbar[i]-1] << std::endl;
					std::cout << pos <<std::endl;
					std::cout << A.getmin() << std::endl;
					std::cout << tjtbar[i] << std::endl;
				 #endif
				 index[tjtbar[i]-1] = pos-1 + A.getmin()[tjtbar[i]-1];
			 } else {
                { 
                int j;
                DimensionIndexIterator<Iterator1D> it(index.getIterator(tjtbar[i],A.getmin()[tjtbar[i]-1],A.getmax()[tjtbar[i]-1]));
				for(j = values.firstIndex(); it.inRange(); it++,j++){
					T eval = abs(A(it.getIndex()) - values(j));
					if(eval > compare){
					  compare = eval;
					  index[tjtbar[i] - 1] = (it.getIndex())[tjtbar[i] - 1];
				    }
			    }
                }
			 }
			 
		  }
	#ifdef _GREEDYDEBUG
		 std::cout << "fatherpivots " << fatherpivots <<  std::endl;
	#endif
		if(fatherpivots.length() > 0){
			values = this->evaluate(index,t,fatherpivots,tjtbarcomplement);
            {
            int j;
			DimensionIndexIterator<IteratorList> it(index.getIterator(tjtbarcomplement,fatherpivots));
            for(j = values.firstIndex(); it.inRange(); it++,++j){
				T eval = abs(A(it.getIndex()) - values(j));
				if(eval > compare){
					compare = eval;
					index.setValue(it.getIndex(),tjtbarcomplement);
				}
			}
            }
		}
	#ifdef _GREEDYDEBUG
		std::cout << "index = " << index << std::endl;
	#endif
	}
	error = abs(A(index) - (*this)(index));
	
	#ifdef _GREEDYDEBUG
	  std::cout << "A(index) = " << A(index) << std::endl;
 	  std::cout << "(*this)(index) = " << (*this)(index) << std::endl;
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, return =  " << index << ", error = " << error << std::endl;
	#endif
	
	return index;
};

template <typename T, typename TensorFunction>
void 
ROA<T, TensorFunction>::approximate(const int rank, const DimensionIndexList &fatherpivots, const int l){
	double error;
	//std::cout << "roa.tcc ROA<T, SVD>::approximate, d = " << d << std::endl;
	DimensionIndex newpiv;
	for(int i = 1; i<=rank; ++i){
		newpiv = GreedyPivotSearch(l,fatherpivots,error);
		if(abs(error) > 0){
		#ifdef _GREEDYDEBUG
				std::cout << "added pivot" << std::endl;
		#endif		
				std::cout << "pivots = " << pivots << std::endl;
			addpivot(newpiv);
		}
	}
};


template <typename T, typename TensorFunction>
void 
ROA<T, TensorFunction>::approximate(const double epsilon, const DimensionIndexList & fatherpivots, const int l){
	double error;
	DimensionIndex greedyres;
	
	do{
		greedyres = GreedyPivotSearch(l,fatherpivots,error);
		#ifdef _GREEDYDEBUG
		  std::cout << "greedyres = " << greedyres << std::endl;
		#endif
		if(pivots.length() > 0){
			if(abs(error/A(pivots[0])) > epsilon){
			   #ifdef _GREEDYDEBUG
				std::cout << "added pivot" << std::endl;
			   #endif			
				addpivot(greedyres);
			}
		} else {
		   #ifdef _GREEDYDEBUG
			std::cout << "added pivot " << std::endl;
		   #endif		
			addpivot(greedyres);
			if(abs(A(pivots[0])) < epsilon){
				break;
			}
		}
	   #ifdef _GREEDYDEBUG
		std::cout << "error/A(pivots[0]) = " << error/A(pivots[0]) << "  epsilon = " << epsilon << ", " << A(pivots[0]) << std::endl;
	   #endif
	} while(abs(error/A(pivots[0])) > epsilon);
};

template <typename T, typename TensorFunction>
void 
ROA<T, TensorFunction>::print() const{
	std::cout << "ROA: d = " << d << "  t = " << t << "   tbar = " << tbar << std::endl;
};

template <typename T, typename TensorFunction>
DimensionIndex 
ROA<T, TensorFunction>::getTbar() const{
return tbar;
};

template <typename T, typename TensorFunction>
DimensionIndex 
ROA<T, TensorFunction>::getT() const{
	return t;
};

template <typename T, typename TensorFunction>
int
ROA<T, TensorFunction>::dim() const{
	return d;
};

template <typename T, typename TensorFunction>
ROA<T,TensorFunction> &
ROA<T,TensorFunction>::operator=(const ROA<T,TensorFunction> &rhs){
	d = rhs.dim();
	t = rhs.getT();
	tbar = rhs.getTbar();
	DimensionIndex tjtbar = t.join(tbar);
	tcomplement = t.getComplement(d);
	A = rhs.A;
	Atensor = rhs.Atensor;
	U = rhs.U;
	V = rhs.V;
	s = rhs.s;
	inverse_s = rhs.inverse_s;
	pivots = rhs.pivots;
	SVD_epsilon = rhs.SVD_epsilon;

	return *this;
}

} // namespace htucker

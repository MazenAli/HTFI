namespace lawa{

template <typename T>
ROA<T, INVERSE>::ROA():d(0),Atensor(1,1),Ainverse(1,1){};

template <typename T>
ROA<T, INVERSE>::ROA(tensor<T> _A, DimensionIndex _t,DimensionIndex _tbar, int _dim):A(_A),t(_t),tbar(_tbar),d(_dim),Atensor(1,1),Ainverse(1,1),tcomplement(_t.getComplement(_dim)){};


template <typename T>
void 
ROA<T, INVERSE>::addpivot(DimensionIndex pivot){
	DimensionIndex ausw;
	this->pivots.add(pivot);
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(pivots.length(),pivots.length());
	if(pivots.length() > 1){
		mat(_(1,pivots.length()-1),_(1,pivots.length()-1)) = Atensor;
	}
	
	Atensor = mat;
	assert(pivot.length() == d);
	ausw = pivots[pivots.length() -1];
	
	for(int i = 0; i < pivots.length(); ++i){
		ausw.setValue(pivots[i],t);
		Atensor(i,pivots.length()) = A(ausw);
	}

	ausw = pivots[pivots.length() -1];

	for(int j = 0; j < pivots.length()-1; ++j){
		ausw.setValue(pivots[i],tcomplement);
		Atensor(pivots.length(),j) = A(ausw);
	}

	Ainverse = Inverse(Atensor);
};

template <typename T>
T 
ROA<T, INVERSE>::operator()(DimensionIndex vals){
	if(pivots.length() == 0) return 0;
	assert(vals.length() == d);
	T sum = 0;

	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > sol(pivots.length(),1); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
	DimensionIndex idx = vals;

	for(int i = 0; i < pivots.length(); ++i){
		idx.setValue(pivots[i],t);
		sol(i,1) = A(idx);
	}
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > refsol;
	refsol = Ainverse*sol;
	sol = refsol;
	
	idx = vals;
	for(int i = 0; i < sol.numRows(); ++i){
		idx.setValue(pivots[i],tcomplement);
		sum+= A(idx)*sol(i + 1,1);
	}
	return sum;
};

template <typename T>
DimensionIndex  
ROA<T, INVERSE>::GreedyPivotSearch(int lmax, DimensionIndexList & fatherpivots, double &error){
	
	  DimensionIndex tjtbar = t.join(tbar);
	  DimensionIndex tjtbarcomplement = tjtbar.getComplement(d);
	  T compare;
	
	  //initialize index with random numbers
	  int fpindex;
	  if(fatherpivots.length()>0){
		  fpindex= ((rand()) % (fatherpivots.length()));
	  }

	  DimensionIndex index(d);
	  index.setRandom(tjtbar,A.minval,A.maxval);
	  index.setRandom(fatherpivots,tjtbarcomplement);
	 
	  compare = abs(A(index) - (*this)(index));

	  
	  for(int l = 1; l <= lmax; ++l){
		  for(int i = 0; i < tjtbar.length(); ++i){
			  for(DimensionIndexIterator<Iterator1D> it = index.getIterator(tjtbar[i],A.minval[i],A.maxval[i]); it.inRange(); it++){
				  T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				  if(eval > compare){
					  compare = eval;
					  index[tjtbar[i] - 1] = (it.getIndex())[tjtbar[i] - 1];
				  }
			  }
		  }

		if(fatherpivots.length() > 0){
			for(DimensionIndexIterator<IteratorList> it = index.getIterator(tjtbarcomplement,fatherpivots); it.inRange(); it++){
				T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				if(eval > compare){
					compare = eval;
					index.setValue(it.getIndex(),tjtbarcomplement);
				}
			}
		}
	}
	error = abs(A(index) - (*this)(index));

	return index;
};

template <typename T>
void 
ROA<T, INVERSE>::approximate(int rank,DimensionIndexList &fatherpivots){
	double error;
	for(int i = 1; i<=rank; ++i){
		addpivot(GreedyPivotSearch(3,fatherpivots,error));
	}
};

template <typename T>
void 
ROA<T, INVERSE>::print(){
	std::cout << "ROA: d = " << d << "  t = " << t << "   tbar = " << tbar << std::endl;
};

template <typename T>
void 
ROA<T, INVERSE>::approximate(double epsilon,DimensionIndexList & fatherpivots){
	double error;
	DimensionIndex greedyres;
	
	do{
		greedyres = GreedyPivotSearch(3,fatherpivots,error);
		if(error > epsilon){
			addpivot(greedyres);
		}
	} while(error > epsilon);
};



//______________________________________________________________________________________________________________________________

template <typename T>
ROA<T, QR>::ROA():d(0),Atensor(1,1),QR(1,1){};

template <typename T>
ROA<T, QR>::ROA(tensor<T> _A, DimensionIndex _t,DimensionIndex _tbar, int _dim):A(_A),t(_t),tbar(_tbar),d(_dim),Atensor(1,1),QR(1,1),tcomplement(_t.getComplement(_dim)){};


template <typename T>
void 
ROA<T, QR>::addpivot(DimensionIndex pivot){
	DimensionIndex ausw;
	this->pivots.add(pivot);
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(pivots.length(),pivots.length());
	if(pivots.length() > 1){
		mat(_(1,pivots.length()-1),_(1,pivots.length()-1)) = Atensor;
	}
	
	Atensor = mat;
	assert(pivot.length() == dim);
	ausw =pivots[pivots.length() -1];
	
	for(int i = 0; i < pivots.length(); ++i){
		ausw.setValue(pivots[i],t);
		Atensor(i,pivots.length()) = A(ausw);
	}

	ausw = pivots[pivots.length() -1];

	for(int j = 0; j < pivots.length()-1; ++j){
		ausw.setValue(pivots[i],tcomplement);
		Atensor(pivots.length(),j) = A(ausw);
	}

	QR = Atensor;
	flens::DenseVector<flens::Array<T> > save;
	qrf(QR,save);
	tauqr = save;
};

template <typename T>
T 
ROA<T, QR>::operator()(DimensionIndex vals){
	if(pivots.length() == 0) return 0;
	assert(vals.length() == d);
	T sum = 0;

	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > sol(pivots.length(),1); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
	DimensionIndex idx = vals;

	for(int i = 0; i < pivots.length(); ++i){
		idx.setValue(pivots[i],t);
		sol(i,1) = A(idx);
	}
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > savesol = sol;
	ormqr(cxxblas::Left,cxxblas::Trans,QR,tauqr,sol);
	trsm(cxxblas::Left,cxxblas::NoTrans,1.,QR.upper(),sol);
	
	idx = vals;
	for(int i = 0; i < sol.numRows(); ++i){
		idx.setValue(pivots[i],tcomplement);
		sum+= A(idx)*sol(i + 1,1);
	}
	return sum;
};

template <typename T>
DimensionIndex  
ROA<T, QR>::GreedyPivotSearch(int lmax, DimensionIndexList & fatherpivots, double &error){
	
	  DimensionIndex tjtbar = t.join(tbar);
	  DimensionIndex tjtbarcomplement = tjtbar.getComplement(d);
	  T compare;
	
	  //initialize index with random numbers
	  int fpindex;
	  if(fatherpivots.length()>0){
		  fpindex= ((rand()) % (fatherpivots.length()));
	  }

	  DimensionIndex index(d);
	  index.setRandom(tjtbar,A.minval,A.maxval);
	  index.setRandom(fatherpivots,tjtbarcomplement);
	 
	  compare = abs(A(index) - (*this)(index));

	  
	  for(int l = 1; l <= lmax; ++l){
		  for(int i = 0; i < tjtbar.length(); ++i){
			  for(DimensionIndexIterator<Iterator1D> it = index.getIterator(tjtbar[i],A.minval[i],A.maxval[i]); it.inRange(); it++){
				  T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				  if(eval > compare){
					  compare = eval;
					  index[tjtbar[i] - 1] = (it.getIndex())[tjtbar[i] - 1];
				  }
			  }
		  }

		if(fatherpivots.length() > 0){
			for(DimensionIndexIterator<IteratorList> it = index.getIterator(tjtbarcomplement,fatherpivots); it.inRange(); it++){
				T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				if(eval > compare){
					compare = eval;
					index.setValue(it.getIndex(),tjtbarcomplement);
				}
			}
		}
	}
	error = abs(A(index) - (*this)(index));

	return index;
};

template <typename T>
void 
ROA<T, QR>::approximate(int rank,DimensionIndexList &fatherpivots){
	double error;
	for(int i = 1; i<=rank; ++i){
		addpivot(GreedyPivotSearch(3,fatherpivots,error));
	}
};

template <typename T>
void 
ROA<T, QR>::print(){
	std::cout << "ROA: d = " << d << "  t = " << t << "   tbar = " << tbar << std::endl;
};

template <typename T>
void 
ROA<T, QR>::approximate(double epsilon,DimensionIndexList & fatherpivots){
	double error;
	DimensionIndex greedyres;
	
	do{
		greedyres = GreedyPivotSearch(3,fatherpivots,error);
		if(error > epsilon){
			addpivot(greedyres);
		}
	} while(error > epsilon);
};


//___________________________________________________________________________________________________________________________________________

template <typename T>
ROA<T, SVD>::ROA():d(0),Atensor(1,1),U(1,1),V(1,1),SVD_epsilon(0.000000000000001){};

template <typename T>
ROA<T, SVD>::ROA(tensor<T> _A, DimensionIndex _t,DimensionIndex _tbar, int _dim):A(_A),t(_t),tbar(_tbar),d(_dim),Atensor(1,1),U(1,1),V(1,1),tcomplement(_t.getComplement(_dim)),SVD_epsilon(0.000000000000001){};


template <typename T>
void 
ROA<T, SVD>::addpivot(DimensionIndex pivot){
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
	svd(Atensave,s,U,V);
	//std::cout << "roa.tcc ROA<T, SVD>::addpivot(DimensionIndex pivot), svd: s = " << s << " U = " << U << " VT = " << V << std::endl;
};

template <typename T>
T 
ROA<T, SVD>::operator()(DimensionIndex vals){
	double eps = SVD_epsilon;
	if(pivots.length() == 0) return 0;
	assert(vals.length() == d);
	T sum = 0;
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > sol(pivots.length(),1); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
	DimensionIndex idx = vals;
	
	for(int i = 0; i < pivots.length(); ++i){
		idx.setValue(pivots[i],t);
		sol(i+1,1) = A(idx);
	}
	
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newdiag(s.length(),s.length());
	for(int i = s.firstIndex(); i<= s.lastIndex(); ++i){
		if(abs(s(i)) > eps){
			newdiag(i,i) = 1/s(i);
		} else {
			newdiag(i,i) = 0;
		}
	}
	
	//std::cout << "roa.tcc operator() Atensor = " << Atensor << std::endl;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > savesol;
	savesol = transpose(U)*sol;
	//flens::blas::mm(cxxblas::Trans,cxxblas::NoTrans,T(1),U,sol,T(0),savesol);
	sol = newdiag * savesol;
	savesol = transpose(V)*sol;
	sol = savesol;
	idx = vals;
	for(int i = 0; i < sol.numRows(); ++i){
		idx.setValue(pivots[i],tcomplement);
		//std::cout << idx << " A(idx) = " << A(idx) << " sol(i + 1,1) = " << sol(i+1,1) << "  " << tcomplement <<   std::endl;
		sum+= A(idx)*sol(i + 1,1);
	}
	return sum;
};

template <typename T>
flens::DenseVector<flens::Array<T> >
ROA<T, SVD>::evaluate(DimensionIndex initialIndex, DimensionIndex _t, int dim, int _min, int _max){
	double eps = SVD_epsilon;
	
	if(pivots.length() == 0){
		flens::DenseVector<flens::Array<T> > ret(_max - _min + 1);
		 for(DimensionIndexIterator<Iterator1D> it = initialIndex.getIterator(dim,_min,_max),int i = 1; it.inRange(); it++,i++){
			ret(i) = A(it.getIndex());
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
		res1 = transpose(U) * Aright;
		for(int i = 1; i <= res1.length();++i){
			if(abs(s(i)) > eps){
				res1(i) = res1(i)/s(i);
			} else {
				res1(i) = 0;
			}
		}

		flens::DenseVector<flens::Array<T> > res2(V.numRows());
		res2 = transpose(V) * res1;
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
		res1 = transpose(Aright) * res2;
	
		return res1;
	}
	
};

template <typename T>
flens::DenseVector<flens::Array<T> >
ROA<T, SVD>::evaluate(DimensionIndex initialIndex, DimensionIndex _t, DimensionIndexList &list, DimensionIndex & list_dims){
	// we assume that all Dimensions in list_dims (father_pivot_dimensions) are not in _t (i.e. they are column indices!)
	double eps = SVD_epsilon;
	if(pivots.length() == 0){
		flens::DenseVector<flens::Array<T> > ret(list.length());
		for(DimensionIndexIterator<IteratorList> it = initialIndex.getIterator(list_dims,list),int i = 1; it.inRange(); it++, i++){
			ret(i) = A(it.getIndex());
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
	res1 = transpose(Aright) * res2;

	return res1;
};

template <typename T>
DimensionIndex  
ROA<T, SVD>::GreedyPivotSearch(int lmax, DimensionIndexList & fatherpivots, double &error){
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, d = " << this->d << std::endl;
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, t = " << t << std::endl;
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tbar = " << tbar << std::endl;
	  DimensionIndex tjtbar = t.join(tbar);
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tjtbar = " << tjtbar << std::endl;
	  DimensionIndex tjtbarcomplement;
	  if(tjtbar.length() < d){
		 tjtbarcomplement = tjtbar.getComplement(d);
	  }
	  T compare;
	
	  //initialize index with random numbers
	  int fpindex;
	  if(fatherpivots.length()>0){
		  fpindex= ((rand()) % (fatherpivots.length()));
	  }

	 
	  DimensionIndex index(d);
	  index.setRandom(tjtbar,A.minval,A.maxval);
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tjbarcomplement.length() = " << tjtbarcomplement.length() << std::endl;
	  if(tjtbarcomplement.length() > 0){
		//index.setRandom(fatherpivots,tjtbarcomplement);
		  index.setValue(fatherpivots[pivots.length() % fatherpivots.length()],tjtbarcomplement);
	  }
	 
	  //std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, index = " << index << std::endl;
	  compare = abs(A(index) - (*this)(index));
	  // std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, compare = " << compare << std::endl;

	  
	  for(int l = 1; l <= lmax; ++l){
		  for(int i = 0; i < tjtbar.length(); ++i){
			  for(DimensionIndexIterator<Iterator1D> it = index.getIterator(tjtbar[i],A.minval[i],A.maxval[i]); it.inRange(); it++){
				  T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				  if(eval > compare){
					  compare = eval;
					  index[tjtbar[i] - 1] = (it.getIndex())[tjtbar[i] - 1];
				  }
			  }
		  }

		if(fatherpivots.length() > 0){
			for(DimensionIndexIterator<IteratorList> it = index.getIterator(tjtbarcomplement,fatherpivots); it.inRange(); it++){
				T eval = abs(A(it.getIndex()) - (*this)(it.getIndex()));
				if(eval > compare){
					compare = eval;
					index.setValue(it.getIndex(),tjtbarcomplement);
				}
			}
		}
	}
	error = abs(A(index) - (*this)(index));
	//std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, return =  " << index << ", error = " << error << std::endl;
	return index;
};

template <typename T>
DimensionIndex  
ROA<T, SVD>::GreedyPivotSearch2(int lmax, DimensionIndexList & fatherpivots, double &error){
	#ifdef _GREEDYDEBUG
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, d = " << this->d << std::endl;
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, t = " << t << std::endl;
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tbar = " << tbar << std::endl;
    #endif
	  DimensionIndex tjtbar = t.join(tbar);
	#ifdef _GREEDYDEBUG 
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tjtbar = " << tjtbar << std::endl;
	#endif
	  DimensionIndex tjtbarcomplement;
	  if(tjtbar.length() < d){
		 tjtbarcomplement = tjtbar.getComplement(d);
	  }
	  T compare;
    #ifdef _GREEDYDEBUG
	  std::cout << "fatherpivots = " << fatherpivots << std::endl;
    #endif
	  DimensionIndex tcomplement = t.getComplement(d);	 
	  DimensionIndex index(d);
	  int count = 0;
	 // do{
		//index.setRandom(tjtbar,A.minval,A.maxval);
		//// std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, tjbarcomplement.length() = " << tjtbarcomplement.length() << std::endl;
		//if(tjtbarcomplement.length() > 0){
		//	index.setRandom(fatherpivots,tjtbarcomplement);
		//	//index.setValue(fatherpivots[pivots.length() % fatherpivots.length()],tjtbarcomplement);
		//}
		//count ++;
	 // }while((( index.equals(pivots,t)) || index.equals(pivots,tcomplement)) && count < 1000);
	  
	  index.setRandom(tjtbar,A.minval,A.maxval);
	  if(tjtbarcomplement.length() > 0){
		  index.setRandom(fatherpivots,tjtbarcomplement);
	  }
	  DimensionIndex saveIndex(d);
	#ifdef _GREEDYDEBUG
	  std::cout << "index = " << index << "A.minval = " << A.minval << " A.maxval = " << A.maxval <<  std::endl;	 
	  std::cout << "A(index) = " << A(index) << std::endl;
	  std::cout << "this->operator()(index) = " << this->operator()(index) << std::endl;
	#endif
	  double startmaximizer = abs(A(index) - this->operator()(index));
	  double savemaximizer;

	  for(int i = 1; i<= 100; i++){
		  saveIndex.setRandom(tjtbar,A.minval,A.maxval);
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
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, index = " << index << std::endl;
	#endif
	  compare = abs(A(index) - (*this)(index));
	#ifdef _GREEDYDEBUG	 
	  std::cout << "roa.tcc ROA<T, SVD>::GreedyPivotSearch, compare = " << compare << std::endl;
	#endif
	  flens::DenseVector<flens::Array<T> > values;
	  for(int l = 1; l <= lmax; ++l){
		 for(int i = 0; i < tjtbar.length(); ++i){
			 values = this->evaluate(index,t,tjtbar[i],A.minval[tjtbar[i]-1],A.maxval[tjtbar[i]-1]);
			 for(DimensionIndexIterator<Iterator1D> it = index.getIterator(tjtbar[i],A.minval[tjtbar[i]-1],A.maxval[tjtbar[i]-1]),int j = values.firstIndex(); it.inRange(); it++,j++){
				 T eval = abs(A(it.getIndex()) - values(j));
				  if(eval > compare){
					  compare = eval;
					  index[tjtbar[i] - 1] = (it.getIndex())[tjtbar[i] - 1];
				  }
			  }
		  }
	#ifdef _GREEDYDEBUG
		 std::cout << "fatherpivots " << std::endl;
	#endif
		if(fatherpivots.length() > 0){
			values = this->evaluate(index,t,fatherpivots,tjtbarcomplement);
			for(DimensionIndexIterator<IteratorList> it = index.getIterator(tjtbarcomplement,fatherpivots), int j = values.firstIndex(); it.inRange(); it++,++j){
				T eval = abs(A(it.getIndex()) - values(j));
				if(eval > compare){
					compare = eval;
					index.setValue(it.getIndex(),tjtbarcomplement);
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

template <typename T>
void 
ROA<T, SVD>::approximate(int rank,DimensionIndexList &fatherpivots){
	double error;
	//std::cout << "roa.tcc ROA<T, SVD>::approximate, d = " << d << std::endl;
	DimensionIndex newpiv;
	for(int i = 1; i<=rank; ++i){
		newpiv = GreedyPivotSearch2(5,fatherpivots,error);
		if(abs(error) > 0){
			addpivot(newpiv);
		}
	}
};

template <typename T>
void 
ROA<T, SVD>::print(){
	std::cout << "ROA: d = " << d << "  t = " << t << "   tbar = " << tbar << std::endl;
};

template <typename T>
void 
ROA<T, SVD>::approximate(double epsilon,DimensionIndexList & fatherpivots){
	double error;
	DimensionIndex greedyres;
	
	do{
		greedyres = GreedyPivotSearch2(3,fatherpivots,error);
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
		}
	   #ifdef _GREEDYDEBUG
		std::cout << "error/A(pivots[0]) = " << error/A(pivots[0]) << "  epsilon = " << epsilon << ", " << A(pivots[0]) << std::endl;
	   #endif
	} while(abs(error/A(pivots[0])) > epsilon);
};


template <typename T>
DimensionIndex 
ROA<T, SVD>::getTbar(){
return tbar;
};

template <typename T>
DimensionIndex 
ROA<T, SVD>::getT(){
	return t;
};

} // namespace lawa

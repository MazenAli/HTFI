

//constructors
template <typename T>
ROA<T>::ROA(tensor<T> _A, flens::DenseVector<flens::Array<int> > _t_index,flens::DenseVector<flens::Array<int> > _tbar_index, int _dim): Atensor(1,1),Ainverse(1,1),UorB(1,1),A(_A), t_index(_t_index),tbar_index(_tbar_index),dim(_dim),is_evaluated(false){
	in_t_index = flens::DenseVector<flens::Array<bool> >(dim);
	in_tbar_index = flens::DenseVector<flens::Array<bool> >(dim);
	for(int i = 1; i <= dim ; ++i){
		in_t_index(i) = false;
		in_tbar_index(i) = false;
		for(int j = t_index.firstIndex(); j <= t_index.lastIndex(); ++j){
			if( i == t_index(j)){
				in_t_index(i) = true;
				break;
			}
		}
		for(int j = tbar_index.firstIndex(); j <= tbar_index.lastIndex(); ++j){
			if( i == tbar_index(j)){
				in_tbar_index(i) = true;
				break;
			}
		}
	}
	
};

template <typename T>
void 
ROA<T>::approximate(int rank,DenseVectorList<int> * fatherpivots){
	double error;
	for(int i = 1; i<=rank; ++i){
		//GreedyPivotSearch(int dimension,tensor<T> A, ROA<T> ROAtTF, flens::DenseVector<flens::Array<int> > indexset, int lmax, DenseVectorListe<int> * fatherpivots){
		addpivot(GreedyPivotSearch(dim,A,*this,t_index,tbar_index,3,fatherpivots,error));
	}
};

template <typename T>
void 
ROA<T>::approximate(double epsilon,DenseVectorList<int> * fatherpivots){
	double error;
	flens::DenseVector<flens::Array<int> > greedyres;
	int foundcounter = 0;
	//std::cout << "Roa parameter" <<  t_index << " " << tbar_index << std::endl;
	do{
		//GreedyPivotSearch(int dimension,tensor<T> A, ROA<T> ROAtTF, flens::DenseVector<flens::Array<int> > indexset, int lmax, DenseVectorListe<int> * fatherpivots){
		greedyres = GreedyPivotSearch(dim,A,*this,t_index,tbar_index,3,fatherpivots,error);
		std::cout << "in ROA<T>::approximate: " << greedyres << " " << error << " funktionswert: " << (A(greedyres) - (*this)(greedyres)) << " A(erg) = " << A(greedyres) << "  this(erg) = " << (*this)(greedyres) <<  std::endl;
		if(error > epsilon){

			//check if this pivot was already found:
			bool found = false;
			if(pivots.length() > 0){
				for(int i = 1; pivots.length() >= i; ++i){
					bool subfound = true;
					for(int j = 1; j<= greedyres.length();++j){
						if( greedyres(j) != (*pivots(i))(j)){
							subfound = false;
							break;
						}
					}
					if(subfound){
						found = true;
						break;
					}
				}
			}
			if(found){
				foundcounter += 1;
			} else {
				std::cout << "added " << std::endl;
				addpivot(greedyres);
			}
		}
		if(foundcounter > 10){
			break;
		}
	} while(error > epsilon*(1+2*sqrt(pivots.length()+0.0)));
};

template <typename T>
void 
ROA<T>::addpivot(flens::DenseVector<flens::Array<int> > pivot){
	flens::DenseVector<flens::Array<int> > ausw(in_t_index.length());
	this->pivots.add(pivot);
	//std::cout << "Atensor resize to " << pivots.length() << std::endl;
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > mat(pivots.length(),pivots.length());
	if(pivots.length() > 1){
		mat(_(1,pivots.length()-1),_(1,pivots.length()-1)) = Atensor;
	}
	//Atensor.resize(pivots.length(),pivots.length(),1,1);
	//std::cout << pivot.length() << "  " << dim << std::endl;
	Atensor = mat;
	assert(pivot.length() == dim);
	for(int i = 1; i <= pivots.length(); ++i){
		for(int k = 1; k <= in_t_index.length(); ++k){
				if(in_t_index(k)){
					ausw(k) = (*pivots(i))(k);
				} else {
					ausw(k) = (*pivots(pivots.length()))(k);
				}
			}
			Atensor(i,pivots.length()) = A(ausw);
		
	}
	for(int j = 1; j <= pivots.length()-1; ++j){
		for(int k = 1; k <= in_t_index.length(); ++k){
				if(in_t_index(k)){
					ausw(k) = (*pivots(pivots.length()))(k);
				} else {
					ausw(k) = (*pivots(j))(k);
				}
			}
			Atensor(pivots.length(),j) = A(ausw);
		
	}
	//

	//boundary methode ist numerisch weniger stabil als FLENS inverse!


	//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > invsave;
	//if(Atensor.numRows() == 1 && Atensor.numCols() == 1){
	//	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > invsave1(1,1);
	//	invsave1(1,1) = 1/Atensor(1,1);
	//	invsave = invsave1;
	//} else {
	//	T alpha = Atensor(Atensor.numRows(),Atensor.numCols());
	//	flens::DenseVector<flens::Array<T> > v;
	//	v = Atensor(_(1,Atensor.numRows()-1),Atensor.numCols());
	//	flens::DenseVector<flens::Array<T> > u;
	//	//std::cout << "here1" << std::endl;
	//	u = Atensor(Atensor.numRows(),_(1,Atensor.numCols()-1));
	//	flens::DenseVector<flens::Array<T> > vbar = Ainverse*v;
	//	flens::DenseVector<flens::Array<T> > ubar = transpose(Ainverse)*u;
	//	
	//	double scalar = alpha - u*vbar;
	//	
	//	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > newinv(Atensor.numRows(),Atensor.numCols());
	//	//std::cout << "vor schleife: " << std::endl;
	//	for(int i = newinv.firstRow(); i<= newinv.lastRow(); ++i){
	//		for(int j = newinv.firstCol(); j <= newinv.lastCol(); ++j){
	//			//std::cout << "i= " << i << " j = " << j << std::endl;
	//			if(i< newinv.lastRow() && j < newinv.lastCol()){
	//				newinv(i,j) = Ainverse(i,j) + (vbar(i)*ubar(j))/scalar;
	//			} else if (i< newinv.lastRow() && j == newinv.lastCol()){
	//				newinv(i,j) = -vbar(i)/scalar;
	//			} else if (i ==  newinv.lastRow() && j < newinv.lastCol()){
	//				newinv(i,j) = -ubar(j)/scalar;
	//			} else {
	//				newinv(i,j) = 1/scalar;
	//			}
	//		}
	//	}
	//	invsave = newinv;
	//}
	//std::cout << "new inverse: " << invsave << std::endl;
	//Ainverse = invsave;
	
	//Ainverse ist auch nicht so gut, zumindest sobald die Anzahl der Pivots wächst..
	//Ainverse = Inverse(Atensor);

	//Versuchen wir die QR-Zerlegung zu berechnen und damit zu lösen..
	
	Ainverse = Atensor;
	flens::DenseVector<flens::Array<T> > save;
	qrf(Ainverse,save);
	tau = save;


	//std::cout << Ainverse << std::endl;
	//std::cout << "Det : " << Det(Atensor) << std::endl;
	//std::cout << "Atensor " << Atensor << std::endl;
	//std::cout << "old inverse: " << Ainverse << std::endl;
	
};

//template <typename T>
//T 
//ROA<T>::operator()(flens::DenseVector<flens::Array<int> > vals){
//	if(pivots.length() == 0) return 0;
//	assert(vals.length() == dim);
//	T sum = 0;
//	//std::cout << "val = " << vals << std::endl;
//	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > sol(pivots.length(),1); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
//	flens::DenseVector<flens::Array<int> > idx(in_t_index.length());
//	//flens::DenseVector<flens::Array<T> > solution;
//	for(int i = 1; i <= pivots.length(); ++i){
//		for(int j = 1; j <= in_t_index.length(); ++j){
//			if(in_t_index(j)){
//				idx(j) = (*pivots(i))(j);	
//			} else {
//				idx(j) = vals(j);
//			}
//		}
//		//std::cout << "idx = " << idx << std::endl;
//	
//		sol(i,1) = A(idx);
//	}
//	//std::cout << "sol = " << sol << std::endl;
//	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > savesol = sol;
//	
//	ormqr(cxxblas::Left,cxxblas::Trans,Ainverse,tau,sol);
//	trsm(cxxblas::Left,cxxblas::NoTrans,1.,Ainverse.upper(),sol);
//
//	//flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > refsol;
//	//refsol = Atensor*sol;
//	//std::cout << "in ROA<T>() vals: " << vals << std::endl;
//	//std::cout << "savesol = " << savesol << std::endl;
//	//std::cout << "refsol = " << refsol << std::endl;
//	
//	for(int i = 1; i <= sol.numRows(); ++i){
//		for(int j = 1; j <= in_t_index.length(); ++j){
//			if(in_t_index(j)){
//				idx(j) = vals(j);	
//			} else {
//				idx(j) = (*pivots(i))(j);
//			}
//		}
//		
//		sum+= A(idx)*sol(i,1);
//	}
//	return sum;
//};

template <typename T>
T 
ROA<T>::operator()(flens::DenseVector<flens::Array<int> > vals){
	if(pivots.length() == 0) return 0;
	assert(vals.length() == dim);
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > res = RankOneApprox(A,*this,vals,1,vals(1),vals(1));
	return res(1,2);
};

template <typename T>
void 
ROA<T>::setUorB(int type,ROA<T> &leftchild,ROA<T> &rightchild){
	if(type == 3){
		//here we should have only one dimension in the t-index, otherwise it is not a leaf.
		assert(t_index.length() == 1);
		UorB.resize(A.maxval(t_index(1)) - A.minval(t_index(1))+1,pivots.length() ,1,1);
		flens::DenseVector<flens::Array<int> > idx(in_t_index.length());
		for(int i = 1; i <= A.maxval(t_index(1)) - A.minval(t_index(1))+1; i++){
			for(int j = 1; j <= pivots.length(); ++j){
				for(int k = 1; k <= in_t_index.length(); ++k){
					if(in_t_index(k)){
						idx(k) = A.minval(t_index(1)) + i - 1;
					} else {
						idx(k) = (*pivots(j))(k);
					}
				}
				//std::cout << idx << std::endl;
				UorB(i,j) = (*this)(idx);
			}
		}
	} else if (type == 1){
		assert(t_index.length() > 1);
	
		flens::DenseVector<flens::Array<T> > tensorprod;
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;
		flens::DenseVector<flens::Array<int> > idx(in_t_index.length());

		UorB.resize(leftchild.pivots.length()*rightchild.pivots.length(),1,1,1);
		for(int i = 1; i <= 1; ++i){
			flens::DenseVector<flens::Array<T> > input(leftchild.pivots.length()*rightchild.pivots.length());
			for(int j = 1; j <= leftchild.pivots.length(); ++j){
				for(int k = 1; k <= rightchild.pivots.length(); ++k){
					for(int l = 1; l <= in_t_index.length(); ++l){
						if(rightchild.in_t_index(l)){
							idx(l) = (*rightchild.pivots(k))(l);
						} else if (leftchild.in_t_index(l)){
						 idx(l) = (*leftchild.pivots(j))(l);
						} else {
							idx(l) = (*pivots(i))(l);
						}
				
					}
					//std::cout << leftchild.in_t_index << std::endl;
					//std::cout << rightchild.in_t_index << std::endl;
					//std::cout << idx << std::endl;
					input((j-1)*rightchild.pivots.length() + k) = A(idx);
				}
			}
			//tensorres = rightchild. * Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			tensorres = Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			ormqr(cxxblas::Left,cxxblas::Trans,rightchild.Ainverse,rightchild.tau,tensorres);
			trsm(cxxblas::Left,cxxblas::NoTrans,1.,rightchild.Ainverse.upper(),tensorres);
			tensorres2 = transpose(tensorres);
			ormqr(cxxblas::Left,cxxblas::Trans,leftchild.Ainverse,leftchild.tau,tensorres2);
			trsm(cxxblas::Left,cxxblas::NoTrans,1.0,leftchild.Ainverse.upper(),tensorres2);
			tensorres = transpose(tensorres2);
			//tensorres2 = tensorres * transpose(leftchild.Ainverse);
			tensorprod = Mat2Vec(tensorres);
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorB(j,i) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
		}

		
	} else {
		assert(t_index.length() > 1);
	
		flens::DenseVector<flens::Array<T> > tensorprod;
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres;
		flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > tensorres2;
		
		flens::DenseVector<flens::Array<int> > idx(in_t_index.length());

		UorB.resize(leftchild.pivots.length()*rightchild.pivots.length(),pivots.length(),1,1);
		for(int i = 1; i <= pivots.length(); ++i){
			flens::DenseVector<flens::Array<T> > input(leftchild.pivots.length()*rightchild.pivots.length());
			for(int j = 1; j <= leftchild.pivots.length(); ++j){
				for(int k = 1; k <= rightchild.pivots.length(); ++k){
					for(int l = 1; l <= in_t_index.length(); ++l){
						if(rightchild.in_t_index(l)){
							idx(l) = (*rightchild.pivots(k))(l);
						} else if (leftchild.in_t_index(l)){
						 idx(l) = (*leftchild.pivots(j))(l);
						} else {
							idx(l) = (*pivots(i))(l);
						}
				
					}
					//std::cout << leftchild.in_t_index << std::endl;
					//std::cout << rightchild.in_t_index << std::endl;
					//std::cout << idx << std::endl;
					input((j-1)*rightchild.pivots.length() + k) = (*this)(idx);
				}
			}
			//std::cout << "childinverse: " << rightchild.Ainverse << std::endl;
			//std::cout << "spalte von A in B: " << input << std::endl;
			//std::cout << "vec2mat : " <<  Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length()) << std::endl;
			tensorres = Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			ormqr(cxxblas::Left,cxxblas::Trans,rightchild.Ainverse,rightchild.tau,tensorres);
			trsm(cxxblas::Left,cxxblas::NoTrans,1.,rightchild.Ainverse.upper(),tensorres);
			//tensorres = rightchild.Ainverse * Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			tensorres2 = transpose(tensorres);
			ormqr(cxxblas::Left,cxxblas::Trans,leftchild.Ainverse,leftchild.tau,tensorres2);
			trsm(cxxblas::Left,cxxblas::NoTrans,1.0,leftchild.Ainverse.upper(),tensorres2);
			tensorres = transpose(tensorres2);
			//tensorres2 = tensorres * transpose(leftchild.Ainverse);
			//std::cout << "tensorres1: " << tensorres2 << std::endl;
			//std::cout << "tensorres2: " << tensorres2 << std::endl;
			tensorprod = Mat2Vec(tensorres);
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorB(j,i) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
		}

	}

};

template <typename T> 
void ROA<T>::evaluateUorB(int type,ROA<T> &leftchild,ROA<T> &rightchild, flens::DenseVector<flens::Array<int> > index){
	if(type == 3){
		//std::cout << "3: " << UorB.numCols() << std::endl;
		evaluate = flens::DenseVector<flens::Array<T> >(UorB.numCols());
		int val = 0;
		for(int i = in_t_index.firstIndex(); i <= in_t_index.lastIndex(); ++i){
			if(in_t_index(i)){
				val = index(i);
			}
		}
		for(int i = 1; i <= UorB.numCols(); ++i){
			evaluate(i) = UorB(val,i);
		}
		is_evaluated = true;
	} else if(type == 1) {
		evaluate = flens::DenseVector<flens::Array<T> >(1);
		//std::cout << "Here1" << std::endl;
		flens::DenseVector<flens::Array <T> > bpart(leftchild.UorB.numCols()*rightchild.UorB.numCols());
		flens::DenseVector<flens::Array <T> > intermediateresult;
		//std::cout << UorB.numRows() << " " << UorB.numCols() << std::endl;
		for(int i= 1; i<= 1; ++i){
			for(int j = 1; j<=leftchild.UorB.numCols()*rightchild.UorB.numCols() ; j++){
				bpart(j) = UorB(j,i);
			}
			//std::cout << "Here2" << std::endl;
			intermediateresult = Vec2Mat(bpart,rightchild.UorB.numCols(),leftchild.UorB.numCols())*leftchild.evaluate;
			//std::cout << "Here4" << std::endl;
			evaluate(i) = 0;
			for(int j = 1; j<= rightchild.UorB.numCols(); ++j){
				evaluate(i) += intermediateresult(j)*rightchild.evaluate(j);
			}
			
		}
		is_evaluated = true;
	} else {
		evaluate = flens::DenseVector<flens::Array<T> >(UorB.numCols());
		flens::DenseVector<flens::Array <T> > bpart(leftchild.UorB.numCols()*rightchild.UorB.numCols());
		flens::DenseVector<flens::Array <T> > intermediateresult;
		
		for(int i= 1; i<= UorB.numCols(); ++i){
			for(int j = 1; j<=leftchild.UorB.numCols()*rightchild.UorB.numCols() ; j++){
				bpart(j) = UorB(j,i);
			}
			
			intermediateresult = Vec2Mat(bpart,rightchild.UorB.numCols(),leftchild.UorB.numCols())*leftchild.evaluate;
			evaluate(i) = 0;
			for(int j = 1; j<= rightchild.UorB.numCols(); ++j){
				evaluate(i) += intermediateresult(j)*rightchild.evaluate(j);
			}
		}
		is_evaluated = true;
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& Stream, ROA<T>& roa)
{
	std::ostream& St = Stream;
	 St << "# of Pivots: " << roa.pivots.length() << std::endl;
	 for(int i = 1; i <= roa.pivots.length(); ++i){
		 St << i << ": " << *(roa.pivots(i)) << std::endl;
	 }
	return St;
};

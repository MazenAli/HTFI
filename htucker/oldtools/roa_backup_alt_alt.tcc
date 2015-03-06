

//constructors
template <typename T>
ROA<T>::ROA(tensor<T> _A, DenseVector<Array<int> > _t_index,DenseVector<Array<int> > _tbar_index, int _dim): Atensor(1,1),Ainverse(1,1),UorB(1,1),A(_A), t_index(_t_index),tbar_index(_tbar_index),dim(_dim),is_evaluated(false){
	in_t_index = DenseVector<Array<bool> >(dim);
	in_tbar_index = DenseVector<Array<bool> >(dim);
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
		//GreedyPivotSearch(int dimension,tensor<T> A, ROA<T> ROAtTF, DenseVector<Array<int> > indexset, int lmax, DenseVectorListe<int> * fatherpivots){
		addpivot(GreedyPivotSearch(dim,A,*this,t_index,tbar_index,3,fatherpivots,error));
	}
};

template <typename T>
void 
ROA<T>::approximate(double epsilon,DenseVectorList<int> * fatherpivots){
	double error;
	DenseVector<Array<int> > greedyres;
	int foundcounter = 0;
	//cout << "Roa parameter" <<  t_index << " " << tbar_index << endl;
	do{
		//GreedyPivotSearch(int dimension,tensor<T> A, ROA<T> ROAtTF, DenseVector<Array<int> > indexset, int lmax, DenseVectorListe<int> * fatherpivots){
		greedyres = GreedyPivotSearch(dim,A,*this,t_index,tbar_index,3,fatherpivots,error);
		cout << "in ROA<T>::approximate: " << greedyres << " " << error << " funktionswert: " << (A(greedyres) - (*this)(greedyres)) << " A(erg) = " << A(greedyres) << "  this(erg) = " << (*this)(greedyres) <<  endl;
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
				cout << "added " << endl;
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
ROA<T>::addpivot(DenseVector<Array<int> > pivot){
	DenseVector<Array<int> > ausw(in_t_index.length());
	this->pivots.add(pivot);
	//cout << "Atensor resize to " << pivots.length() << endl;
	GeMatrix<FullStorage<T,ColMajor> > mat(pivots.length(),pivots.length());
	if(pivots.length() > 1){
		mat(_(1,pivots.length()-1),_(1,pivots.length()-1)) = Atensor;
	}
	//Atensor.resize(pivots.length(),pivots.length(),1,1);
	//cout << pivot.length() << "  " << dim << endl;
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


	//GeMatrix<FullStorage<T,ColMajor> > invsave;
	//if(Atensor.numRows() == 1 && Atensor.numCols() == 1){
	//	GeMatrix<FullStorage<T,ColMajor> > invsave1(1,1);
	//	invsave1(1,1) = 1/Atensor(1,1);
	//	invsave = invsave1;
	//} else {
	//	T alpha = Atensor(Atensor.numRows(),Atensor.numCols());
	//	DenseVector<Array<T> > v;
	//	v = Atensor(_(1,Atensor.numRows()-1),Atensor.numCols());
	//	DenseVector<Array<T> > u;
	//	//cout << "here1" << endl;
	//	u = Atensor(Atensor.numRows(),_(1,Atensor.numCols()-1));
	//	DenseVector<Array<T> > vbar = Ainverse*v;
	//	DenseVector<Array<T> > ubar = transpose(Ainverse)*u;
	//	
	//	double scalar = alpha - u*vbar;
	//	
	//	GeMatrix<FullStorage<T,ColMajor> > newinv(Atensor.numRows(),Atensor.numCols());
	//	//cout << "vor schleife: " << endl;
	//	for(int i = newinv.firstRow(); i<= newinv.lastRow(); ++i){
	//		for(int j = newinv.firstCol(); j <= newinv.lastCol(); ++j){
	//			//cout << "i= " << i << " j = " << j << endl;
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
	//cout << "new inverse: " << invsave << endl;
	//Ainverse = invsave;
	
	//Ainverse ist auch nicht so gut, zumindest sobald die Anzahl der Pivots wächst..
	//Ainverse = Inverse(Atensor);

	//Versuchen wir die QR-Zerlegung zu berechnen und damit zu lösen..
	
	Ainverse = Atensor;
	DenseVector<Array<T> > save;
	qrf(Ainverse,save);
	tau = save;


	//cout << Ainverse << endl;
	//cout << "Det : " << Det(Atensor) << endl;
	//cout << "Atensor " << Atensor << endl;
	//cout << "old inverse: " << Ainverse << endl;
	
};

//template <typename T>
//T 
//ROA<T>::operator()(DenseVector<Array<int> > vals){
//	if(pivots.length() == 0) return 0;
//	assert(vals.length() == dim);
//	T sum = 0;
//	//cout << "val = " << vals << endl;
//	GeMatrix<FullStorage<T,ColMajor> > sol(pivots.length(),1); //first hat die spaltenindizes aus dem pivot, last die zeilenindizes
//	DenseVector<Array<int> > idx(in_t_index.length());
//	//DenseVector<Array<T> > solution;
//	for(int i = 1; i <= pivots.length(); ++i){
//		for(int j = 1; j <= in_t_index.length(); ++j){
//			if(in_t_index(j)){
//				idx(j) = (*pivots(i))(j);	
//			} else {
//				idx(j) = vals(j);
//			}
//		}
//		//cout << "idx = " << idx << endl;
//	
//		sol(i,1) = A(idx);
//	}
//	//cout << "sol = " << sol << endl;
//	GeMatrix<FullStorage<T,ColMajor> > savesol = sol;
//	
//	ormqr(Left,Trans,Ainverse,tau,sol);
//	trsm(Left,NoTrans,1.,Ainverse.upper(),sol);
//
//	//GeMatrix<FullStorage<T,ColMajor> > refsol;
//	//refsol = Atensor*sol;
//	//cout << "in ROA<T>() vals: " << vals << endl;
//	//cout << "savesol = " << savesol << endl;
//	//cout << "refsol = " << refsol << endl;
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
ROA<T>::operator()(DenseVector<Array<int> > vals){
	if(pivots.length() == 0) return 0;
	assert(vals.length() == dim);
	GeMatrix<FullStorage<T,ColMajor> > res = RankOneApprox(A,*this,vals,1,vals(1),vals(1));
	return res(1,2);
};

template <typename T>
void 
ROA<T>::setUorB(int type,ROA<T> &leftchild,ROA<T> &rightchild){
	if(type == 3){
		//here we should have only one dimension in the t-index, otherwise it is not a leaf.
		assert(t_index.length() == 1);
		UorB.resize(A.maxval(t_index(1)) - A.minval(t_index(1))+1,pivots.length() ,1,1);
		DenseVector<Array<int> > idx(in_t_index.length());
		for(int i = 1; i <= A.maxval(t_index(1)) - A.minval(t_index(1))+1; i++){
			for(int j = 1; j <= pivots.length(); ++j){
				for(int k = 1; k <= in_t_index.length(); ++k){
					if(in_t_index(k)){
						idx(k) = A.minval(t_index(1)) + i - 1;
					} else {
						idx(k) = (*pivots(j))(k);
					}
				}
				//cout << idx << endl;
				UorB(i,j) = (*this)(idx);
			}
		}
	} else if (type == 1){
		assert(t_index.length() > 1);
	
		DenseVector<Array<T> > tensorprod;
		GeMatrix<FullStorage<T,ColMajor> > tensorres;
		GeMatrix<FullStorage<T,ColMajor> > tensorres2;
		DenseVector<Array<int> > idx(in_t_index.length());

		UorB.resize(leftchild.pivots.length()*rightchild.pivots.length(),1,1,1);
		for(int i = 1; i <= 1; ++i){
			DenseVector<Array<T> > input(leftchild.pivots.length()*rightchild.pivots.length());
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
					//cout << leftchild.in_t_index << endl;
					//cout << rightchild.in_t_index << endl;
					//cout << idx << endl;
					input((j-1)*rightchild.pivots.length() + k) = A(idx);
				}
			}
			//tensorres = rightchild. * Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			tensorres = Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			ormqr(Left,Trans,rightchild.Ainverse,rightchild.tau,tensorres);
			trsm(Left,NoTrans,1.,rightchild.Ainverse.upper(),tensorres);
			tensorres2 = transpose(tensorres);
			ormqr(Left,Trans,leftchild.Ainverse,leftchild.tau,tensorres2);
			trsm(Left,NoTrans,1.0,leftchild.Ainverse.upper(),tensorres2);
			tensorres = transpose(tensorres2);
			//tensorres2 = tensorres * transpose(leftchild.Ainverse);
			tensorprod = Mat2Vec(tensorres);
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorB(j,i) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
		}

		
	} else {
		assert(t_index.length() > 1);
	
		DenseVector<Array<T> > tensorprod;
		GeMatrix<FullStorage<T,ColMajor> > tensorres;
		GeMatrix<FullStorage<T,ColMajor> > tensorres2;
		
		DenseVector<Array<int> > idx(in_t_index.length());

		UorB.resize(leftchild.pivots.length()*rightchild.pivots.length(),pivots.length(),1,1);
		for(int i = 1; i <= pivots.length(); ++i){
			DenseVector<Array<T> > input(leftchild.pivots.length()*rightchild.pivots.length());
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
					//cout << leftchild.in_t_index << endl;
					//cout << rightchild.in_t_index << endl;
					//cout << idx << endl;
					input((j-1)*rightchild.pivots.length() + k) = (*this)(idx);
				}
			}
			//cout << "childinverse: " << rightchild.Ainverse << endl;
			//cout << "spalte von A in B: " << input << endl;
			//cout << "vec2mat : " <<  Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length()) << endl;
			tensorres = Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			ormqr(Left,Trans,rightchild.Ainverse,rightchild.tau,tensorres);
			trsm(Left,NoTrans,1.,rightchild.Ainverse.upper(),tensorres);
			//tensorres = rightchild.Ainverse * Vec2Mat(input,rightchild.pivots.length(),leftchild.pivots.length());
			tensorres2 = transpose(tensorres);
			ormqr(Left,Trans,leftchild.Ainverse,leftchild.tau,tensorres2);
			trsm(Left,NoTrans,1.0,leftchild.Ainverse.upper(),tensorres2);
			tensorres = transpose(tensorres2);
			//tensorres2 = tensorres * transpose(leftchild.Ainverse);
			//cout << "tensorres1: " << tensorres2 << endl;
			//cout << "tensorres2: " << tensorres2 << endl;
			tensorprod = Mat2Vec(tensorres);
			for(int j = tensorprod.firstIndex(); j <= tensorprod.lastIndex(); ++j){
				UorB(j,i) = tensorprod(j);  // Spaltenweise eintragen... fertig!
			}
		}

	}

};

template <typename T> 
void ROA<T>::evaluateUorB(int type,ROA<T> &leftchild,ROA<T> &rightchild, DenseVector<Array<int> > index){
	if(type == 3){
		//cout << "3: " << UorB.numCols() << endl;
		evaluate = DenseVector<Array<T> >(UorB.numCols());
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
		evaluate = DenseVector<Array<T> >(1);
		//cout << "Here1" << endl;
		DenseVector<Array <T> > bpart(leftchild.UorB.numCols()*rightchild.UorB.numCols());
		DenseVector<Array <T> > intermediateresult;
		//cout << UorB.numRows() << " " << UorB.numCols() << endl;
		for(int i= 1; i<= 1; ++i){
			for(int j = 1; j<=leftchild.UorB.numCols()*rightchild.UorB.numCols() ; j++){
				bpart(j) = UorB(j,i);
			}
			//cout << "Here2" << endl;
			intermediateresult = Vec2Mat(bpart,rightchild.UorB.numCols(),leftchild.UorB.numCols())*leftchild.evaluate;
			//cout << "Here4" << endl;
			evaluate(i) = 0;
			for(int j = 1; j<= rightchild.UorB.numCols(); ++j){
				evaluate(i) += intermediateresult(j)*rightchild.evaluate(j);
			}
			
		}
		is_evaluated = true;
	} else {
		evaluate = DenseVector<Array<T> >(UorB.numCols());
		DenseVector<Array <T> > bpart(leftchild.UorB.numCols()*rightchild.UorB.numCols());
		DenseVector<Array <T> > intermediateresult;
		
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
ostream& operator<<(ostream& Stream, ROA<T>& roa)
{
	ostream& St = Stream;
	 St << "# of Pivots: " << roa.pivots.length() << endl;
	 for(int i = 1; i <= roa.pivots.length(); ++i){
		 St << i << ": " << *(roa.pivots(i)) << endl;
	 }
	return St;
};


namespace lawa{

kscdomatrix::kscdomatrix(int _j0, int _J, double _alpha, double _beta, double _Tstep, DenseVector<Array<double> > _intensities):alpha(_alpha), beta(_beta), Tstep(_Tstep),j0(_j0), J(_J),intensities(_intensities),cumWavNumber(_J-_j0+1){
	transitionintensities ti(intensities);
	double eps = 0.000000001;
	tensor<double> lambdatens(ti,intensities.length(),0,3);
	HTuckerTree<double,SVD> lambdatree(intensities.length(),lambdatens);
	lambdatree.approximate(eps);
	setUorB(lambdatree);

	emptytensor et;
	tensor<double> empty(et,intensities.length(),0,1);
	HTuckerTree<double,SVD> ones(intensities.length(),empty);

	for(GeneralTreeIterator<HTuckerTreeNode<double,SVD> > TIT1 = ones.getGeneralTree().end(); TIT1 >= ones.getGeneralTree().begin();TIT1--){
		if(TIT1.getNode()->isLeaf()){
			GeMatrix<FullStorage<double,ColMajor> > tmp(2,1);
			tmp = 1,
				  1;
			TIT1.getNode()->getContent()->setUorB(tmp);
			TIT1.getNode()->getContent()->UorB_numel = 1;
		} else {
			GeMatrix<FullStorage<double,ColMajor> > tmp(1,1);
			tmp = 1.0;
			TIT1.getNode()->getContent()->setUorB(tmp);
			TIT1.getNode()->getContent()->UorB_numel = 1;
			TIT1.getNode()->getContent()->UorB_rcnumel = 1;
			TIT1.getNode()->getContent()->UorB_lcnumel = 1;
		}
	}

	HTuckerTree<double,SVD> lambdacum = lambdatree*ones;

	transitionmatrix tm(lambdatree,lambdacum);
	tensor<double> transmat(tm,intensities.length(),0,3);
	lambda = HTuckerTree<double,SVD>(intensities.length(),transmat);
	lambda.approximate(eps);
	setUorB(lambda);


	Basis<double,Orthogonal,Interval,Multi> mwbasis(3,j0);
	mwbasis.enforceBoundaryCondition<DirichletBC>();

	int N = mwbasis.mra.cardI(J);

	DimensionIndex maxmw = factor(N);
	DimensionIndex minmw(1,maxmw.length());
	DimensionIndex minlambda(0,intensities.length());
	DimensionIndex maxlambda(1,intensities.length());

	maxvals = maxlambda.join(maxmw);
	minvals = minlambda.join(minmw);
	
	cout << "N = " << N << endl;
	cout << "maxmw = "  << maxmw << endl;
	cout << " maxvals " << maxvals << endl;
	
	for(int i = j0; i <= J; ++i){
		cumWavNumber[i-j0] = mwbasis.mra.cardI(i);
	}

	cout << "cumWavNumber = " << cumWavNumber << endl;


}

double 
kscdomatrix::operator() (DimensionIndex vals){
	DimensionIndex wav(vals.length() - intensities.length());
	DimensionIndex ints(intensities.length());
	for(int i = 0; i < vals.length() - intensities.length(); ++i){
		wav[i] = vals[i + intensities.length()];
	}
	for(int i = 0; i < intensities.length(); ++i){
		ints[i] = vals[i];
	}
	 
	// Zerlege Index in Zeilen - und Spaltenindex
	DimensionIndex zeilewav(wav.length());
	DimensionIndex spaltewav(wav.length());

	for(int i = 0; i < wav.length(); ++i){
		spaltewav[i] = floor((wav[i]-minvals[intensities.length()+i]+0.0)/(0.0 + maxvals[intensities.length()+i]-minvals[intensities.length()+i]+1))+minvals[intensities.length()+i];
		zeilewav[i] = wav[i]-minvals[intensities.length()+i] - (spaltewav[i] - minvals[intensities.length()])*(maxvals[intensities.length()]-minvals[intensities.length()]+1);
	}
	
	cout << spaltewav << endl;
	cout << zeilewav << endl;
	
	int zeilewert = 0;
	int spaltewert = 0;
	for(int i = 0; i < vals.length() - intensities.length() ; ++i){
		if(i > 0){
			zeilewert *= maxvals[intensities.length() + i] - minvals[intensities.length() + i] + 1;
			spaltewert *= maxvals[intensities.length() + i] - minvals[intensities.length() + i] + 1;
		}
		zeilewert = zeilewert + (zeilewav[i] - minvals[intensities.length() + i]);
		spaltewert = spaltewert + (spaltewav[i] - minvals[intensities.length() + i]);
	
	}
	zeilewert = zeilewert + 1;
	spaltewert += 1;
	
	int zeilelevel = 0;
	int spaltelevel = 0;

	for(int i = 0; i < cumWavNumber.length(); ++i){
		if(i == 0){
			if(zeilewert > 0 && zeilewert <= cumWavNumber[i]){
				zeilelevel = i;
			}
			if(spaltewert > 0 && spaltewert <= cumWavNumber[i]){
				spaltelevel = i;
			}
		} else {
			if(zeilewert > cumWavNumber[i-1] && zeilewert <= cumWavNumber[i]){
				zeilelevel = i;
			}
			if(spaltewert > cumWavNumber[i-1] && spaltewert <= cumWavNumber[i]){
				spaltelevel = i;
			}
		}
	}
	int zeilenverschiebung = 0;
	int spaltenverscheibung = 0;

	if(zeilelevel > 0){
		zeilenverschiebung = zeilewert - cumWavNumber[zeilelevel-1];
	} else {
		zeilenverschiebung = zeilewert;
	}

	if(spaltelevel > 0){
		spaltenverscheibung = spaltewert - cumWavNumber[spaltelevel-1];
	} else {
		spaltenverscheibung = spaltewert;
	}


	if(zeilelevel == 0 && spaltelevel == 0){
	
	} else if (zeilelevel == 0) {
	
	} else if (spaltelevel == 0){
	
	} else {
	
	}
	return 1.0;
}



transitionintensities::transitionintensities(DenseVector<Array<double> > _intensities): intensities(_intensities),_min(DimensionIndex(_intensities.length())),_max(DimensionIndex(_intensities.length())),_inner(DimensionIndex(_intensities.length())),zeilendims(DimensionIndex(_intensities.length())),spaltendims(DimensionIndex(_intensities.length())),zeilendimscum(DimensionIndex(_intensities.length())),spaltendimscum(DimensionIndex(_intensities.length())){
	int prod = 1;
	//cout << "here " << endl;
	for(int i = 1; i <= intensities.length(); ++i){
		//cerr << " i = " << i << endl;
		_min[i-1] = 0;
		_max[i-1] = 3;
		_inner[i-1] = sqrt(_max[i-1] - _min[i-1] + 1.0)+_min[i-1]-1;
		prod *= (_max[i-1] - _min[i-1])+1;
		spaltendims[i-1] = _inner[i-1] - _min[i-1] + 1;
		zeilendims[i-1] = (_max[i-1] - _min[i-1] + 1)/spaltendims[i-1];
		spaltendimscum[i-1] = 1;
		zeilendimscum[i-1] = 1;
		for(int j = 0; j < i-1; ++j){
			spaltendimscum[j]*=spaltendims[i-1];
			zeilendimscum[j]*=zeilendims[i-1];
		}
	}
}

double
transitionintensities::operator ()(lawa::DimensionIndex vals){
	
	assert(vals.length() == _min.length());
	int zeile,spalte,z;
	
	for(int i = 0; i < _min.length(); ++i){
		if(i==0){
			zeile = ((vals[i] - _min[i]) % zeilendims[i])*zeilendimscum[i];
			spalte = (((vals[i] - _min[i]) - zeile + 1)/zeilendims[i])*spaltendimscum[i];
		} else {
			z = ((vals[i] - _min[i]) % zeilendims[i]);
			zeile += z*spaltendimscum[i];
			spalte +=  ((int) ((vals[i] - _min[i]) - z)/zeilendims[i]) * zeilendimscum[i];
		}
	}
	if(zeile == spalte){
		return 0;
	} else if ((abs(log((zeile ^ spalte)+0.0)/log(2.0) - floor(log((zeile ^ spalte)+0.0)/log(2.0)))) < 0.00000001 && spalte > zeile) {
		return intensities(floor(log((zeile ^ spalte) +0.0)/log(2.0)+1));
	} else {
		return 0;
	}
}


transitionmatrix::transitionmatrix(HTuckerTree<double,SVD> &_t1, HTuckerTree<double,SVD> &_t2): t1(_t1),t2(_t2),_min(DimensionIndex(_t1.dim())),_max(DimensionIndex(_t1.dim())),_inner(DimensionIndex(_t1.dim())),zeilendims(DimensionIndex(_t1.dim())),spaltendims(DimensionIndex(_t1.dim())),zeilendimscum(DimensionIndex(_t1.dim())),spaltendimscum(DimensionIndex(_t1.dim())){
	int prod = 1;
	//cout << "here " << endl;
	cout << t1.dim() << endl;
	t1.print_w_UorB();
	for(int i = 1; i <= t1.dim(); ++i){
		//cerr << " i = " << i << endl;
		_min[i-1] = 0;
		_max[i-1] = 3;
		_inner[i-1] = sqrt(_max[i-1] - _min[i-1] + 1.0)+_min[i-1]-1;
		prod *= (_max[i-1] - _min[i-1])+1;
		spaltendims[i-1] = _inner[i-1] - _min[i-1] + 1;
		zeilendims[i-1] = (_max[i-1] - _min[i-1] + 1)/spaltendims[i-1];
		spaltendimscum[i-1] = 1;
		zeilendimscum[i-1] = 1;
		for(int j = 0; j < i-1; ++j){
			spaltendimscum[j]*=spaltendims[i-1];
			zeilendimscum[j]*=zeilendims[i-1];
		}
	}
}

double
transitionmatrix::operator ()(lawa::DimensionIndex vals){
	
	assert(vals.length() == _min.length());
	int zeile,spalte,z;
	DimensionIndex zeilenindex(t1.dim());
	
	for(int i = 0; i < _min.length(); ++i){
		if(i==0){
			zeilenindex[i] = ((vals[i] - _min[i]) % zeilendims[i]);
			zeile = ((vals[i] - _min[i]) % zeilendims[i])*zeilendimscum[i];
			spalte = (((vals[i] - _min[i]) - zeile + 1)/zeilendims[i])*spaltendimscum[i];
		} else {
			z = ((vals[i] - _min[i]) % zeilendims[i]);
			zeilenindex[i] = z;
			zeile += z*spaltendimscum[i];
			spalte +=  ((int) ((vals[i] - _min[i]) - z)/zeilendims[i]) * zeilendimscum[i];
		}
	}
	if(zeile == spalte){
		//cout << "t2.evaluate( " << zeilenindex << ")" << endl;
		return t2.evaluate(zeilenindex);
	} else {
		//cout << "t1.evaluate(" << vals << ")" << endl;
		return t1.evaluate(vals);
	}
}


}//namespace lawa
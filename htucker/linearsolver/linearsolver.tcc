namespace htucker{

#include <cmath>

template <typename T>
int
cg(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol, const long maxIterations){
	
	T alpha, beta, rNorm, rNormPrev, bNorm;
	bNorm = b.L2norm();
	double eps = 0.000000000001;
	//std::cout << "tol = " << tol << std::endl;
	//std::cout << " A = ";
	//A.print_values();
	//std::cout << " x = ";
	//x.print_values();
	HTuckerTree<T> r = A*x;
	//std::cout << "r = A*x = ";
	//r.print_values();
	r = reapproximate(r,eps);
	//std::cout << "reapprox(r) = ";
	//r.print_values();
	r = b - r;
	//std::cout << "r = b - r = ";
	//r.print_values();
	r = reapproximate(r,eps);
	//std::cout << "reapprox(r) = ";
	//r.print_values();

	//b.print_values();

    HTuckerTree<T> p = r;
    rNorm = r.L2norm();
    for (long k=1; k<=maxIterations; k++) {
		std::cout << " in htuckercg k = " << k << std::endl;
		std::cout << " in htuckercg normr = " << rNorm << std::endl;
        
        if (rNorm/bNorm<=tol) {
            return k-1;
        }
        HTuckerTree<T> Ap = A*p;
		//std::cout << "A*p = ";
		//Ap.print_values();
		Ap = reapproximate(Ap,eps);
		//std::cout << "reapproximate(Ap) = ";
		//Ap.print_values();
		//Ap.print_w_UorB();
		//std::cout << "p = ";
		//p.print_values();
		//p.print_w_UorB();
		//std::cout << "(p.ScalarProduct(Ap)) = " << (p.ScalarProduct(Ap)) << std::endl;
		//std::cout << "rNorm*rNorm = " << rNorm*rNorm << std::endl;
        alpha = rNorm*rNorm/(p.ScalarProduct(Ap));
        //std::cout << "alpha = " << alpha << std::endl;
		x = x + alpha*p;
		//std::cout << "x = x+alpha * p = ";
		//x.print_values();
		x = reapproximate(x,eps);
        //std::cout << "x = reapproximate(x) = ";
		//x.print_values(); 
		r = r - alpha*Ap;
		//std::cout << " r = r - alpha*Ap = ";
		//r.print_values();
		r = reapproximate(r,eps);
		//r.print_values();
        rNormPrev = rNorm;
        
		rNorm = r.L2norm();
        //std::cout << "rNorm = " << rNorm << std::endl;
		beta = rNorm*rNorm/(rNormPrev*rNormPrev);
       // std::cout << "beta = " << beta << std::endl;
		p = beta*p + r;
		//std::cout << " p = beta*p + r = ";
		//p.print_values();
		p = reapproximate(p,eps);
		//std::cout << " p = reapproximate(p) = ";
		//p.print_values();
    }
    return maxIterations;
	return 1;
}


template <typename T>
int
cg(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol, const long maxIterations){
	
	T alpha, beta, rNorm, rNormPrev, bNorm;
	bNorm = b.L2norm();
	double eps = 0.000000000001;
	//std::cout << "tol = " << tol << std::endl;
	//std::cout << " A = ";
	//A.print_values();
	//std::cout << " x = ";
	//x.print_values();
	HTuckerTree<T> r = A*x;
	//std::cout << "r = A*x = ";
	//r.print_values();
	r = truncate(r,maxrank);
	//std::cout << "reapprox(r) = ";
	//r.print_values();
	r = b - r;
	//std::cout << "r = b - r = ";
	//r.print_values();
	r = truncate_nonorthogonal(r,maxrank);
	//std::cout << "reapprox(r) = ";
	//r.print_values();

	//b.print_values();

    HTuckerTree<T> p = r;
    rNorm = r.L2norm();
    for (long k=1; k<=maxIterations; k++) {
		std::cout << " in htuckercg k = " << k << std::endl;
		std::cout << " in htuckercg normr = " << rNorm << std::endl;
		std::cout << " error " << rNorm/bNorm << std::endl;
        
        if (rNorm/bNorm<=tol) {
            return k-1;
        }
        HTuckerTree<T> Ap = A*p;
		//std::cout << "A*p = ";
		//Ap.print_values();
		Ap = truncate(Ap,maxrank);
		//std::cout << "reapproximate(Ap) = ";
		//Ap.print_values();
		//Ap.print_w_UorB();
		//std::cout << "p = ";
		//p.print_values();
		//p.print_w_UorB();
		//std::cout << "(p.ScalarProduct(Ap)) = " << (p.ScalarProduct(Ap)) << std::endl;
		//std::cout << "rNorm*rNorm = " << rNorm*rNorm << std::endl;
        alpha = rNorm*rNorm/(p.ScalarProduct(Ap));
        //std::cout << "alpha = " << alpha << std::endl;
		x = x + alpha*p;
		//std::cout << "x = x+alpha * p = ";
		//x.print_values();
		x = truncate_nonorthogonal(x,maxrank);
        //std::cout << "x = reapproximate(x) = ";
		//x.print_values(); 
		r = r - alpha*Ap;
		//std::cout << " r = r - alpha*Ap = ";
		//r.print_values();
		r = truncate_nonorthogonal(r,maxrank);
		//r.print_values();
        rNormPrev = rNorm;
        
		rNorm = r.L2norm();
        //std::cout << "rNorm = " << rNorm << std::endl;
		beta = rNorm*rNorm/(rNormPrev*rNormPrev);
       // std::cout << "beta = " << beta << std::endl;
		p = beta*p + r;
		//std::cout << " p = beta*p + r = ";
		//p.print_values();
		p = truncate_nonorthogonal(p,maxrank);
		//std::cout << " p = reapproximate(p) = ";
		//p.print_values();
    }
    return maxIterations;
	return 1;
}





template <typename T>
int
bicgstab(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol, const long maxIterations){
	
	T eps = tol/10;
	HTuckerTree<T> r = A*x;
	
	r = reapproximate(r,eps);
	r = b - r;
	r = reapproximate(r,eps);
	
	//b.print_values();
	T rho = 1;
	T alpha = 1;
	T omega = 1;


	HTuckerTree<T> p;
	HTuckerTree<T> r0 = r;
	HTuckerTree<T> v;

	if(abs(sqrt(r0.L2norm())) < tol){
		std::cerr << "initial guess is too close to 0" << std::endl;
		return 1;
	}

    p = r;
    
    for (long k=1; k<=maxIterations; k++) {
		T save_rho = rho;
		rho = r0.ScalarProduct(r);
		T beta = rho/save_rho*alpha/omega;
		if(k == 1){
			p = r0;
		} else {
			HTuckerTree<T> temp = p-omega*v;
			temp = reapproximate(temp, eps);
			p = r + beta*temp;
			p = reapproximate(p,eps);
		}
		v = A*p;
		v = reapproximate(v, eps);
		alpha = rho/r0.ScalarProduct(v);
		HTuckerTree<T> s = r - alpha*v;
		s= reapproximate(s, eps);
		HTuckerTree<T> t = A*s;
		t = reapproximate(t, eps);
		omega = t.ScalarProduct(s)/t.ScalarProduct(t);
		HTuckerTree<T> temp2 = alpha*p + omega*s;
		temp2 = reapproximate(temp2, eps);
		x = x + temp2;
		x = reapproximate(x, eps);
		if(std::fabs(alpha)+std::fabs(omega) < eps){
			break;
		}

		r = s-omega*t;
		r = reapproximate(r, eps);
    }
    return maxIterations;
	return 1;
}






template <typename T, typename Aop, class A1, class A2>
int
bicgstab(const HTuckerClosure<Aop,A1,A2> &A,  HTuckerTree<T> &x,const HTuckerTree<T> &b, const T tol, const long maxIterations){

	//std::cout << "Hallo Welt " << std::endl;
    T eps = tol/10;
	
	HTuckerTree<T> r = A*x;

	//std::cout << " r =>>>>>>>>>>>>>>>>" << std::endl;
	//r.print_w_UorB(); 
	r = reapproximate(r,eps);
	r = b - r;
	r = reapproximate(r,eps);
	
	//b.print_values();
	T rho = 1;
	T alpha = 1;
	T omega = 1;


	HTuckerTree<T> p;
	HTuckerTree<T> r0 = r;
	HTuckerTree<T> v;

	if(abs(sqrt(r0.L2norm())) < tol){
		std::cerr << "initial guess is too close to 0" << std::endl;
		return 1;
	}

    p = r;
    
    for (long k=1; k<=maxIterations; k++) {
		T save_rho = rho;
		rho = r0.ScalarProduct(r);
		T beta = rho/save_rho*alpha/omega;
		if(k == 1){
			p = r0;
		} else {
			HTuckerTree<T> temp = p-omega*v;
			temp = reapproximate(temp, eps);
			p = r + beta*temp;
			p = reapproximate(p,eps);
		}
		//p.print_w_UorB();
		v = A*p;

		v = reapproximate(v, eps);
		alpha = rho/r0.ScalarProduct(v);
		HTuckerTree<T> s = r - alpha*v;
		s= reapproximate(s, eps);
		HTuckerTree<T> t = A*s;
		t = reapproximate(t, eps);
		omega = t.ScalarProduct(s)/t.ScalarProduct(t);
		HTuckerTree<T> temp2 = alpha*p + omega*s;
		temp2 = reapproximate(temp2, eps);
		HTuckerTree<T> temp3 = x + temp2;
		x = reapproximate(temp3, eps);
		
		r = s-omega*t;
		r = reapproximate(r, eps);
			
		T error = r.L2norm();
		#ifdef SOLVER_DEBUG
			//std::cout << "iteration: " << k << std::endl;
			//std::cout << "fehler 2 = " << error << std::endl;
		#endif
		if(error < eps){
			#ifdef SOLVER_DEBUG
std:std::cout << "termination after " << k << "iterations error = " << error << std::endl;
			#endif
			break;
		}
		
    }
    return maxIterations;
	return 1;
}

template <typename T>
int
bicgstab(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol, const long maxIterations){
	
	T eps = tol/10;
	HTuckerTree<T> r = A*x;
	
	r = truncate(r,maxrank);
	r = b - r;
	r = truncate_nonorthogonal(r,maxrank);
	
	//b.print_values();
	T rho = 1;
	T alpha = 1;
	T omega = 1;


	HTuckerTree<T> p;
	HTuckerTree<T> r0 = r;
	HTuckerTree<T> v;

	if(abs(sqrt(r0.L2norm())) < tol){
		std::cerr << "initial guess is too close to 0" << std::endl;
		return 1;
	}

    p = r;
    
    for (long k=1; k<=maxIterations; k++) {
		T save_rho = rho;
		rho = r0.ScalarProduct(r);
		T beta = rho/save_rho*alpha/omega;
		if(k == 1){
			p = r0;
		} else {
			HTuckerTree<T> temp = p-omega*v;
			temp = truncate_nonorthogonal(temp, maxrank);
			p = r + beta*temp;
			p = truncate_nonorthogonal(p,maxrank);
		}
		v = A*p;
		//v.print_values();
		v = truncate(v, maxrank);
		//v.print_values();
		alpha = rho/r0.ScalarProduct(v);
		HTuckerTree<T> s = r - alpha*v;
		//s.print_values();
		
		s= truncate_nonorthogonal(s, maxrank);
		//s.print_values();
		HTuckerTree<T> t = A*s;
		t = truncate(t, maxrank);
		omega = t.ScalarProduct(s)/t.ScalarProduct(t);
		HTuckerTree<T> temp2 = alpha*p + omega*s;
		//temp2.print_w_UorB();
		//temp2.print_values();
		temp2 = truncate_nonorthogonal(temp2, maxrank);
		//temp2.print_w_UorB();
		temp2.print_values();
		x = x + temp2;
		x = truncate_nonorthogonal(x, maxrank);
		//std::cout << "error = " << temp2.L2norm() << std::endl;
		if(temp2.L2norm() < eps){
			break;
		}

		r = s-omega*t;
		r = truncate_nonorthogonal(r, maxrank);
    }
    return maxIterations;
	return 1;
}

template <typename T, typename Aop, class A1, class A2>
int
bicgstab(const HTuckerClosure<Aop,A1,A2> &A,  HTuckerTree<T> &x,const HTuckerTree<T> &b, const int maxrank, const T tol, const long maxIterations){

	//std::cout << "Hallo Welt " << std::endl;
    T eps = tol/10;
	
	HTuckerTree<T> r = A*x;

	//std::cout << " r =>>>>>>>>>>>>>>>>" << std::endl;
	//r.print_w_UorB(); 
	//r.print_values();
	r = truncate_nonorthogonal(r,maxrank);
	//r.print_values();
	r = b - r;
	//r.print_values();
	r = truncate_nonorthogonal(r,maxrank);
	//r.orthogonalize();
	//r.print_values();
	//std::cout << "hier 1 " << std::endl;
	//b.print_values();
	T rho = 1;
	T alpha = 1;
	T omega = 1;


	HTuckerTree<T> p;
	HTuckerTree<T> r0 = r;
	HTuckerTree<T> v;

	if(abs(sqrt(r0.L2norm())) < tol){
		std::cerr << "initial guess is too close to 0" << std::endl;
		return 1;
	}

    p = r;
    
    for (long k=1; k<=maxIterations; k++) {
		T save_rho = rho;
		rho = r0.ScalarProduct(r);
		T beta = rho/save_rho*alpha/omega;
		if(k == 1){
			p = r0;
		} else {
			HTuckerTree<T> temp = p-omega*v;
			//temp.print_values();
			temp = truncate_nonorthogonal(temp, maxrank);
			//temp.print_values();
			p = r + beta*temp;
			//p.print_values();
			p = truncate_nonorthogonal(p,maxrank);
			//p.orthogonalize();
			//p.print_values();
		}
		//std::cout << "hier 2 " << std::endl;
		//p.print_w_UorB();
		v = A*p;
		//v.print_values();
		v = truncate(v, maxrank);
		//v.print_values();
		//std::cout << "hier 3 " << std::endl;
		alpha = rho/r0.ScalarProduct(v);
		HTuckerTree<T> s = r - alpha*v;
		//s.print_values();
		s= truncate_nonorthogonal(s, maxrank);
		//s.orthogonalize();
		//s.print_values();
		HTuckerTree<T> t = A*s;
		//std::cout << "hier 4 " << std::endl;
		//t.print_values();
		t = truncate_nonorthogonal(t, maxrank);
		//t.print_values();
		omega = t.ScalarProduct(s)/t.ScalarProduct(t);
		HTuckerTree<T> temp2 = alpha*p + omega*s;
		//temp2.print_values();
		temp2 = truncate_nonorthogonal(temp2, maxrank);
		//temp2.orthogonalize();
		//temp2.print_values();
		//std::cout << "hier 5 " << std::endl;
		HTuckerTree<T> temp3 = x + temp2;
		//temp3.print_values();
		x = truncate_nonorthogonal(temp3, maxrank);
		//std::cout << "hier 6" << std::endl;
		//x.print_values();
		r = s-omega*t;
		r = truncate_nonorthogonal(r, maxrank);
		//r.orthogonalize();
		//std::cout << "hier 7" << std::endl;
			
		T error = r.L2norm();
		
		if(error < eps){
			#ifdef SOLVER_DEBUG
			std::cout << "iteration: " << k << std::endl;
			std::cout << "fehler 2 = " << error << std::endl;
			std::cout << "increment = " << temp2.L2norm() << std::endl;
			#endif
			break;
		}
		
    }
    return maxIterations;
	return 1;
}




template <typename T, typename Aop, class A1, class A2>
int
bicgstablinf(const HTuckerClosure<Aop,A1,A2> &A,  HTuckerTree<T> &x,const HTuckerTree<T> &b, const T tol, const long maxIterations){

	//std::cout << "Hallo Welt " << std::endl;
    T eps = tol/10;
	
	HTuckerTree<T> r = A*x;

	//std::cout << " r =>>>>>>>>>>>>>>>>" << std::endl;
	//r.print_w_UorB(); 
	r = reapproximate(r,eps);
	r = b - r;
	r = reapproximate(r,eps);
	
	//b.print_values();
	T rho = 1;
	T alpha = 1;
	T omega = 1;


	HTuckerTree<T> p;
	HTuckerTree<T> r0 = r;
	HTuckerTree<T> v;

	if(abs(sqrt(r0.L2norm())) < tol){
		std::cerr << "initial guess is too close to 0" << std::endl;
		return 1;
	}

    p = r;
    
    for (long k=1; k<=maxIterations; k++) {
		T save_rho = rho;
		rho = r0.ScalarProduct(r);
		T beta = rho/save_rho*alpha/omega;
		if(k == 1){
			p = r0;
		} else {
			HTuckerTree<T> temp = p-omega*v;
			temp = reapproximate(temp, eps);
			p = r + beta*temp;
			p = reapproximate(p,eps);
		}
		//p.print_w_UorB();
		v = A*p;

		v = reapproximate(v, eps);
		alpha = rho/r0.ScalarProduct(v);
		HTuckerTree<T> s = r - alpha*v;
		s= reapproximate(s, eps);
		HTuckerTree<T> t = A*s;
		t = reapproximate(t, eps);
		omega = t.ScalarProduct(s)/t.ScalarProduct(t);
		HTuckerTree<T> temp2 = alpha*p + omega*s;
		temp2 = reapproximate(temp2, eps);
		HTuckerTree<T> temp3 = x + temp2;
		x = reapproximate(temp3, eps);
		
		r = s-omega*t;
		r = reapproximate(r, eps);
			
		T error = r.L2norm();
		#ifdef SOLVER_DEBUG
			std::cout << "iteration: " << k << std::endl;
			std::cout << "fehler 2 = " << error << std::endl;
		#endif
		if(error < eps){
			break;
		}
		
    }
    return maxIterations;
	return 1;
}







template <typename T>
int
gmres(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol, long maxIterations){
    using flens::_;

	double eps = tol/100.0;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	
	#ifdef SOLVER_DEBUG
		std::cerr << " N = " << N << std::endl;
	#endif

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

	HTuckerTree<T> prod = A*x;
	prod = reapproximate(prod,eps);
	HTuckerTree<T> r = b - prod;
	r = reapproximate(r,eps);
	T beta = r.L2norm();
	#ifdef SOLVER_DEBUG
		std::cerr << " beta = " << beta << std::endl;
	#endif

    if (beta==0.0) {
        return 0;
    }

	HTuckerTree<T> *V = (HTuckerTree<T> *) malloc ((maxIterations + 1)*sizeof(HTuckerTree<T>));
	HTuckerTree<T> *tmp = new (V) HTuckerTree<T>(A.dim());
    V[0] = r *(1/ beta);
      
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >  H(maxIterations+1, maxIterations);
   

	HTuckerTree<T> w_j = r;
	flens::DenseVector<flens::Array<T> > g(maxIterations+1);
	flens::DenseVector<flens::Array<T> > c(maxIterations+1), s(maxIterations+1);

	T nu, rho = tol + 1;
	T Htemp, h_ij;
    
	g(1) = beta;
    
    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }


		w_j =  A * V[j-1];
		w_j = reapproximate(w_j,eps);
		T normInitialWj = w_j.L2norm();
		for (int i=1; i<=j; ++i) {
			H(i,j) = w_j.ScalarProduct(V[i-1]);
			w_j = w_j - H(i,j) * V[i-1];
			#ifdef SOLVER_DEBUG
				w_j.print_values();
			#endif // SOLVER_DEBUG
			w_j = reapproximate(w_j,eps);
			#ifdef SOLVER_DEBUG
				w_j.print_values();
				w_j.print_w_UorB();
			#endif // SOLVER_DEBUG
		}
		H(j+1,j) = w_j.L2norm();

		if (H(j+1,j) / normInitialWj < 1.0) {
			for (int i=1; i<=j; ++i) {
				#ifdef SOLVER_DEBUG
					V[i-1].print_w_UorB();
				#endif // SOLVER_DEBUG
          		Htemp = w_j.ScalarProduct(V[i-1]);
		   		w_j = w_j - Htemp * V[i-1];
				w_j = reapproximate(w_j,eps);
			}	
			H(j+1, j) = w_j.L2norm();
		}

		for (int i=1; i<=j-1; ++i) {
			h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
			H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
			H(i,j) =    h_ij;
		}

		nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
		if (nu!=0.0) {
			HTuckerTree<T> *tmp = new (V+j) HTuckerTree<T>(A.dim());
		
			V[j] = w_j *(1/H(j+1,j));
		}

		if (nu!=0.0) {
			s(j) = H(j+1,j)/nu;
			c(j) = H(j,j)/nu;
			H(j,j) = nu;
			H(j+1,j) = 0.0;
			g(j+1) = -s(j)*g(j);
			g(j) = c(j)*g(j);
		}
		rho = fabs(g(j+1)) / beta;
	}

        
    if (j>=1) {
		flens::DenseVector<flens::Array<T> > y(j);
		for (int i=j; i>=1; --i) {
			y(i) = g(i) / H(i,i);
			for (int l=j; l>i; --l) {
				y(i) -= H(i,l) * y(l) / H(i,i);
			}
		}
		
		for(int i = 1; i <= y.length(); ++i){
			x = x +  V[i-1] * y(i);
			x = reapproximate(x,eps);
		}
	}


	free(V);   
	return j;
}

template <typename T>
int
gmres(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol, long maxIterations){
    using flens::_;

	double eps = tol/100.0;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	
	#ifdef SOLVER_DEBUG
		std::cerr << " N = " << N << std::endl;
	#endif

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

	HTuckerTree<T> prod = A*x;
	prod = truncate(prod,maxrank);
	HTuckerTree<T> r = b - prod;
	r = truncate_nonorthogonal(r,maxrank);
	T beta = r.L2norm();
	#ifdef SOLVER_DEBUG
		//std::cerr << " beta = " << beta << std::endl;
	#endif

    if (beta==0.0) {
        return 0;
    }

	HTuckerTree<T> *V = (HTuckerTree<T> *) malloc ((maxIterations + 1)*sizeof(HTuckerTree<T>));
	HTuckerTree<T> *tmp = new (V) HTuckerTree<T>(A.dim());
    V[0] = r *(1/ beta);
      
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >  H(maxIterations+1, maxIterations);
   

	HTuckerTree<T> w_j = r;
	flens::DenseVector<flens::Array<T> > g(maxIterations+1);
	flens::DenseVector<flens::Array<T> > c(maxIterations+1), s(maxIterations+1);

	T nu, rho = tol + 1;
	T Htemp, h_ij;
    
	g(1) = beta;
    
    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }


		w_j =  A * V[j-1];
		w_j = truncate(w_j,maxrank);
		T normInitialWj = w_j.L2norm();
		for (int i=1; i<=j; ++i) {
			H(i,j) = w_j.ScalarProduct(V[i-1]);
			w_j = w_j - H(i,j) * V[i-1];
			#ifdef SOLVER_DEBUG
				//w_j.print_values();
			#endif // SOLVER_DEBUG
			w_j = truncate_nonorthogonal(w_j,maxrank);
			#ifdef SOLVER_DEBUG
				//w_j.print_values();
				//w_j.print_w_UorB();
			#endif // SOLVER_DEBUG
		}
		H(j+1,j) = w_j.L2norm();

		if (H(j+1,j) / normInitialWj < 1.0) {
			for (int i=1; i<=j; ++i) {
				#ifdef SOLVER_DEBUG
					//V[i-1].print_w_UorB();
				#endif // SOLVER_DEBUG
          		Htemp = w_j.ScalarProduct(V[i-1]);
		   		w_j = w_j - Htemp * V[i-1];
				w_j = truncate_nonorthogonal(w_j,maxrank);
			}	
			H(j+1, j) = w_j.L2norm();
		}

		for (int i=1; i<=j-1; ++i) {
			h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
			H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
			H(i,j) =    h_ij;
		}

		nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
		if (nu!=0.0) {
			HTuckerTree<T> *tmp = new (V+j) HTuckerTree<T>(A.dim());
		
			V[j] = w_j *(1/H(j+1,j));
		}

		if (nu!=0.0) {
			s(j) = H(j+1,j)/nu;
			c(j) = H(j,j)/nu;
			H(j,j) = nu;
			H(j+1,j) = 0.0;
			g(j+1) = -s(j)*g(j);
			g(j) = c(j)*g(j);
		}
		rho = fabs(g(j+1)) / beta;
	}

        
    if (j>=1) {
		flens::DenseVector<flens::Array<T> > y(j);
		for (int i=j; i>=1; --i) {
			y(i) = g(i) / H(i,i);
			for (int l=j; l>i; --l) {
				y(i) -= H(i,l) * y(l) / H(i,i);
			}
		}
		
		for(int i = 1; i <= y.length(); ++i){
			x = x +  V[i-1] * y(i);
			x = truncate_nonorthogonal(x,maxrank);
		}
	}


	free(V);   
	return j;
}

template <typename T>
int
gmresm(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol,
       const int restart, long maxIterations)
{
    
	double eps = tol/100.;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
	HTuckerTree<T> prod = A*x;
	prod = reapproximate(prod,eps);
	HTuckerTree<T> r = b-prod;
	r = reapproximate(r,eps);
    while (r.L2norm()>tol && k<=maxIterations) {
        k+=gmres(A,x,b,tol,restart);
        prod = A*x;
		prod = reapproximate(prod,eps);
		r=prod -b;
		r = reapproximate(r,eps);
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << r.L2norm() << std::endl;
        #endif
    }
    return k;
}

template <typename T>
int
gmresm(const HTuckerTree<T> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol,
       const int restart, long maxIterations)
{
    
	double eps = tol/100.;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
	HTuckerTree<T> prod = A*x;
	prod = truncate(prod,maxrank);
	HTuckerTree<T> r = b-prod;
	r = truncate_nonorthogonal(r,maxrank);
    while (r.L2norm()>tol && k<=maxIterations) {
        k+=gmres(A,x,b,maxrank,tol,restart);
        prod = A*x;
		prod = truncate(prod,maxrank);
		r=prod -b;
		r = truncate_nonorthogonal(r,maxrank);
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << r.L2norm() << std::endl;
        #endif
    }
    return k;
}

template <typename T, typename Aop, class A1, class A2>
int
gmres(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol, long maxIterations){
    using flens::_;

	double eps = tol/100.0;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	
	#ifdef SOLVER_DEBUG
		std::cerr << " N = " << N << std::endl;
	#endif

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

	HTuckerTree<T> prod = A*x;
	prod = reapproximate(prod,eps);
	HTuckerTree<T> r = b - prod;
	r = reapproximate(r,eps);
	T beta = r.L2norm();
	#ifdef SOLVER_DEBUG
		std::cerr << " beta = " << beta << std::endl;
	#endif

    if (beta==0.0) {
        return 0;
    }

	HTuckerTree<T> *V = (HTuckerTree<T> *) malloc ((maxIterations + 1)*sizeof(HTuckerTree<T>));
	HTuckerTree<T> *tmp = new (V) HTuckerTree<T>(A.dim());
    V[0] = r *(1/ beta);
      
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >  H(maxIterations+1, maxIterations);
   

	HTuckerTree<T> w_j = r;
	flens::DenseVector<flens::Array<T> > g(maxIterations+1);
	flens::DenseVector<flens::Array<T> > c(maxIterations+1), s(maxIterations+1);

	T nu, rho = tol + 1;
	T Htemp, h_ij;
    
	g(1) = beta;
    
    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }


		w_j =  A * V[j-1];
		w_j = reapproximate(w_j,eps);
		T normInitialWj = w_j.L2norm();
		for (int i=1; i<=j; ++i) {
			H(i,j) = w_j.ScalarProduct(V[i-1]);
			w_j = w_j - H(i,j) * V[i-1];
			#ifdef SOLVER_DEBUG
				w_j.print_values();
			#endif // SOLVER_DEBUG
			w_j = reapproximate(w_j,eps);
			#ifdef SOLVER_DEBUG
				w_j.print_values();
				w_j.print_w_UorB();
			#endif // SOLVER_DEBUG
		}
		H(j+1,j) = w_j.L2norm();

		if (H(j+1,j) / normInitialWj < 1.0) {
			for (int i=1; i<=j; ++i) {
				#ifdef SOLVER_DEBUG
					V[i-1].print_w_UorB();
				#endif // SOLVER_DEBUG
          		Htemp = w_j.ScalarProduct(V[i-1]);
		   		w_j = w_j - Htemp * V[i-1];
				w_j = reapproximate(w_j,eps);
			}	
			H(j+1, j) = w_j.L2norm();
		}

		for (int i=1; i<=j-1; ++i) {
			h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
			H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
			H(i,j) =    h_ij;
		}

		nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
		if (nu!=0.0) {
			HTuckerTree<T> *tmp = new (V+j) HTuckerTree<T>(A.dim());
		
			V[j] = w_j *(1/H(j+1,j));
		}

		if (nu!=0.0) {
			s(j) = H(j+1,j)/nu;
			c(j) = H(j,j)/nu;
			H(j,j) = nu;
			H(j+1,j) = 0.0;
			g(j+1) = -s(j)*g(j);
			g(j) = c(j)*g(j);
		}
		rho = fabs(g(j+1)) / beta;
	}

        
    if (j>=1) {
		flens::DenseVector<flens::Array<T> > y(j);
		for (int i=j; i>=1; --i) {
			y(i) = g(i) / H(i,i);
			for (int l=j; l>i; --l) {
				y(i) -= H(i,l) * y(l) / H(i,i);
			}
		}
		
		for(int i = 1; i <= y.length(); ++i){
			x = x +  V[i-1] * y(i);
			x = reapproximate(x,eps);
		}
	}


	free(V);   
	return j;
}


template <typename T, typename Aop, class A1, class A2>
int
gmres(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol, long maxIterations){
    using flens::_;

	double eps = tol/100.0;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}
	
	#ifdef SOLVER_DEBUG
		std::cerr << " N = " << N << std::endl;
	#endif

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

	HTuckerTree<T> prod = A*x;
	prod = truncate(prod,maxrank);
	HTuckerTree<T> r = b - prod;
	r = truncate_nonorthogonal(r,maxrank);
	T beta = r.L2norm();
	#ifdef SOLVER_DEBUG
		std::cerr << " beta = " << beta << std::endl;
	#endif

    if (beta==0.0) {
        return 0;
    }

	HTuckerTree<T> *V = (HTuckerTree<T> *) malloc ((maxIterations + 1)*sizeof(HTuckerTree<T>));
	HTuckerTree<T> *tmp = new (V) HTuckerTree<T>(A.dim());
    V[0] = r *(1/ beta);
      
	flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> >  H(maxIterations+1, maxIterations);
   

	HTuckerTree<T> w_j = r;
	flens::DenseVector<flens::Array<T> > g(maxIterations+1);
	flens::DenseVector<flens::Array<T> > c(maxIterations+1), s(maxIterations+1);

	T nu, rho = tol + 1;
	T Htemp, h_ij;
    
	g(1) = beta;
    
    int j;
    for (j=0; j<=maxIterations-1;) {
        #ifdef SOLVER_DEBUG
            std::cerr << "j = " << j << ", rho = " << rho << std::endl;
        #endif
        if (rho <= tol) {
            break;
        } else {
               ++j;
        }


		w_j =  A * V[j-1];
		w_j = truncate(w_j,maxrank);
		T normInitialWj = w_j.L2norm();
		for (int i=1; i<=j; ++i) {
			H(i,j) = w_j.ScalarProduct(V[i-1]);
			w_j = w_j - H(i,j) * V[i-1];
			#ifdef SOLVER_DEBUG
				//w_j.print_values();
			#endif // SOLVER_DEBUG
			w_j = truncate_nonorthogonal(w_j,maxrank);
			#ifdef SOLVER_DEBUG
				//w_j.print_values();
				//w_j.print_w_UorB();
			#endif // SOLVER_DEBUG
		}
		H(j+1,j) = w_j.L2norm();

		if (H(j+1,j) / normInitialWj < 1.0) {
			for (int i=1; i<=j; ++i) {
				#ifdef SOLVER_DEBUG
				//V[i-1].print_w_UorB();
				#endif // SOLVER_DEBUG
          		Htemp = w_j.ScalarProduct(V[i-1]);
		   		w_j = w_j - Htemp * V[i-1];
				w_j = truncate_nonorthogonal(w_j,maxrank);
			}	
			H(j+1, j) = w_j.L2norm();
		}

		for (int i=1; i<=j-1; ++i) {
			h_ij =      c(i) * H(i,j) + s(i) * H(i+1,j);
			H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j);
			H(i,j) =    h_ij;
		}

		nu = sqrt(H(_(j,j+1),j) * H(_(j,j+1),j));
		if (nu!=0.0) {
			HTuckerTree<T> *tmp = new (V+j) HTuckerTree<T>(A.dim());
		
			V[j] = w_j *(1/H(j+1,j));
		}

		if (nu!=0.0) {
			s(j) = H(j+1,j)/nu;
			c(j) = H(j,j)/nu;
			H(j,j) = nu;
			H(j+1,j) = 0.0;
			g(j+1) = -s(j)*g(j);
			g(j) = c(j)*g(j);
		}
		rho = fabs(g(j+1)) / beta;
	}

        
    if (j>=1) {
		flens::DenseVector<flens::Array<T> > y(j);
		for (int i=j; i>=1; --i) {
			y(i) = g(i) / H(i,i);
			for (int l=j; l>i; --l) {
				y(i) -= H(i,l) * y(l) / H(i,i);
			}
		}
		
		for(int i = 1; i <= y.length(); ++i){
			x = x +  V[i-1] * y(i);
			x = truncate_nonorthogonal(x,eps);
		}
	}


	free(V);   
	return j;
}

template <typename T, typename Aop, class A1, class A2>
int
gmresm(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const T tol,
       const int restart, long maxIterations)
{
    
	double eps = tol/100.;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
	HTuckerTree<T> prod = A*x;
	prod = reapproximate(prod,eps);
	HTuckerTree<T> r = b-prod;
	r = reapproximate(r,eps);
    while (r.L2norm()>tol && k<=maxIterations) {
        k+=gmres(A,x,b,tol,restart);
        prod = A*x;
		prod = reapproximate(prod,eps);
		r=prod -b;
		r = reapproximate(r,eps);
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << r.L2norm() << std::endl;
        #endif
    }
    return k;
}

template <typename T, typename Aop, class A1, class A2>
int
gmresm(const HTuckerClosure<Aop,A1,A2> &A, HTuckerTree<T> &x, const HTuckerTree<T> &b, const int maxrank, const T tol,
       const int restart, long maxIterations)
{
    
	double eps = tol/100.;
    int N = 1;

	for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = b.getGeneralTree().begin(); TIT <= b.getGeneralTree().end(); TIT++){
		if(TIT.getNode()->isLeaf()){
			N*=TIT.getNode()->getContent()->getUorB().numRows();
		}
	}

    if (maxIterations==-1) {
        maxIterations=2*N;
    }

    int k=0;
	HTuckerTree<T> prod = A*x;
	prod = truncate(prod,maxrank);
	HTuckerTree<T> r = b-prod;
	r = truncate_nonorthogonal(r,maxrank);
    while (r.L2norm()>tol && k<=maxIterations) {
        k+=gmres(A,x,b,maxrank,tol,restart);
        prod = A*x;
		prod = truncate(prod,maxrank);
		r=prod -b;
		r = truncate_nonorthogonal(r,maxrank);
        #ifdef SOLVER_DEBUG
            std::cerr << "k = " << k << ", rho = " << r.L2norm() << std::endl;
        #endif
    }
    return k;
}


} // namespace htucker

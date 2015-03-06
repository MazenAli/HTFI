namespace  htucker{


	
	template <typename T>
	HTuckerTree<T> positivepart(const HTuckerTree<T> & tree,const int maxrank, const T tol, const T newtontol, long maxIterations){
		HTuckerTree<T> result = 0.5*(tree+absolutevalue(tree,maxrank,tol,newtontol,maxIterations));
		result = truncate_nonorthogonal(result,maxrank);
		return result;
	}

	template <typename T>
	HTuckerTree<T> absolutevalue(const HTuckerTree<T> & tree,const int maxrank, const T tol, const T newtontol, long maxIterations){
		HTuckerTree<T> squaretree = vec2diag(tree)*tree;
		//squaretree.print_values();
		squaretree = truncate_nonorthogonal(squaretree,maxrank);
		DimensionIndex idx(tree.dim());
		for(GeneralTreeIterator<HTuckerTreeNode<T> > TIT = squaretree.getGeneralTree().end(); TIT >= squaretree.getGeneralTree().begin(); TIT--){
			if(TIT.getNode()->isLeaf()){
				idx[(TIT.getNode()->getContent()->getIndex())[0]-1] = TIT.getNode()->getContent()->getUorB().numRows();
			}
		}


		HTuckerTree<T> x_n = squaretree*0.5;
		HTuckerTree<T> squarediag = vec2diag(squaretree);
		//squaretree.print_values();
		HTuckerTree<T> one = ones<T>(idx);
		HTuckerTree<T> inverse = one;

		for(long i = 1; i <= maxIterations; ++i){
			HTuckerTree<T> secondterm = vec2diag(x_n);
			inverse = one;
			std::cout << cg(secondterm,inverse,one,maxrank,tol,maxIterations) << std::endl;
			HTuckerTree<T> x_n_old = x_n;
			x_n = squarediag*inverse;
			//x_n.print_w_UorB();
			x_n = truncate_nonorthogonal(x_n,maxrank);
			HTuckerTree<T> error = x_n - x_n_old;
			x_n = 0.5*x_n + 0.5*x_n_old;
			std::cout << i << ": error = " << error.L2norm() << std::endl;
			if(error.L2norm() < newtontol){
				return x_n;
			}
			x_n = truncate_nonorthogonal(x_n,maxrank);
			

		}


		return x_n;
	}

} //namespace htucker

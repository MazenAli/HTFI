#ifndef HTUCKER_LAWAEXT_PLOTFUNCTION_H
#define HTUCKER_LAWAEXT_PLOTFUNCTION_H 1

namespace lawa{

	template <typename _T, FunctionSide _Side, DomainType _Domain, Construction _Cons>
	void
	plotfunction(const BasisFunction<_T,_Side,_Domain,_Cons> * BF, const int j, const long k, const unsigned short deriv, int steps = 100){
		flens::DenseVector<flens::Array< _T> > y;
		flens::DenseVector<flens::Array< _T> > x = linspace(BF->support(j,k).l1,BF->support(j,k).l2,steps+1);
		lawa::evaluatefunction(BF,j,k,deriv,steps,y);
		
		Gnuplot::set_GNUPlotPath("C:\\Users\\Andi\\Documents\\Programmieren\\gnuplot\\gnuplot\\bin");
		Gnuplot g1 = Gnuplot("pm3d");
		std::cout << x.firstIndex() << "  " << x.lastIndex() << std::endl;
		std::cout << y.firstIndex() << " " << y.lastIndex() << std::endl;
		std::ofstream plotFileT("plotBasisFunction.tmp");
		for(int i = 1; i <= x.length(); ++i){
			plotFileT << x(i-1) << " " << y(i) << std::endl;
		} 
		plotFileT.close();

		g1.cmd("plot 'plotBasisFunction.tmp' with lines");
		g1.cmd("set terminal png");
		g1.cmd("set output 'plotBasisFunction.png' ");
		g1.cmd("replot");
	}


}

#endif // HTUCKER_LAWAEXT_PLOTFUNCTION_H

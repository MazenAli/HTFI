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

#ifndef LAWA_PLOTTOOLS_H
#define LAWA_PLOTTOOLS_H 1

#include <flens/flens.h>
#include <lawa/lawa.h>

using namespace flens;
namespace lawa{

template <typename _T, FunctionSide _Side, DomainType _Domain, Construction _Cons>
void
evaluate(BasisFunction<_T,_Side,_Domain,_Cons> * BF, int j, long k, unsigned short deriv, int steps, DenseVector<Array<_T> > &x){
	assert(steps>=2);
	DenseVector<Array<_T> > y(steps+1);
	Support<_T> supp = BF->support(j,k);
	_T dist = supp.l2 - supp.l1;
	for(int i = 1; i<= steps+1; ++i){
		y(i) = BF->operator()(supp.l1+(i-1)*dist/(steps+0.0),j,k,deriv);
	}
	x = y;
}





} //namespace lawa

#endif // LAWA_PLOTTOOLS_H

/*
  Copyright (c) 2010-12 Andreas Paffenholz
 
  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2, or (at your option) any
  later version: http://www.gnu.org/licenses/gpl.txt.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/PowerSet.h"
#include "polymake/Graph.h"
#include "polymake/Set.h"
#include "polymake/graph/connected.h"

namespace polymake { namespace bimatrix_games {


    PowerSet<int>
    connected_components ( const Array<std::pair<Vector<Rational>,Vector<Rational> > > & nash ) {
     
      Set<std::pair<Set<int>, Set<int> > > cc;

      Graph<> G(nash.size());
      int i = 0;
      Array<std::pair<Vector<Rational>,Vector<Rational> > >::const_iterator nit = nash.begin();
      Array<std::pair<Vector<Rational>,Vector<Rational> > >::const_iterator nit_end = nash.end(); nit_end--;
      for ( ; nit != nit_end; ++nit, ++i ) {
	Array<std::pair<Vector<Rational>,Vector<Rational> > >::const_iterator n2it = nit; ++n2it;
	int j = i+1;
	for ( ; n2it != nash.end(); ++n2it, ++j ) 
	  if (  (*nit).first == (*n2it).first ||  (*nit).second == (*n2it).second )
	    G.edge(i,j);
      }
      
      return graph::connected_components(G);
    }


    UserFunction4perl("# @category Game Theory", 
		      &connected_components, "connected_components(Array<Pair<Vector<Rational>,Vector<Rational>>>)");

  }}

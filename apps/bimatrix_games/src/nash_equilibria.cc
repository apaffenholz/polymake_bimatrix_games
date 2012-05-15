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
#include "polymake/Matrix.h"
#include "polymake/SparseMatrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/Set.h"
#include "polymake/Array.h"

namespace polymake { namespace bimatrix_games {


    /***************************************************

 compute the nash equilibria of a non-generate game 
 with feasible reduced best response polyhedra using the following algorithm:

 Let A, B be the payoff matrices
 we search for solutions of the following system: 
 
 Ay+r=\1
 B^tx+s=\1
 x^tr+y^ts=0
 x,y,r,s\ge 0
 
Compute vertices of the polyhedron p defined by Ay\le \1, y\ge 0
for each vertex v of P compute vertices of the polyhedron q defined by B^tx\le \1, x\ge 0,  x^tr=y^ts=0
for each vertex w of q the pair (v,w) is a nash equilibrium

    ******************************************************/

    // FIXME: choose role of A and B depending on row/column dimensions
    Array<std::pair<Vector<Rational>,Vector<Rational> > >
    nash_equilibria ( const std::pair<Matrix<Rational>,Matrix<Rational> > & PM ) {
     
      Set<std::pair<Vector<Rational>,Vector<Rational> > > nash;   // collects the extreme nash equilibria

      perl::Object p("polytope::Polytope<Rational>");

      Matrix<Rational> A = ones_vector<Rational>(PM.first.rows())|-PM.first;
      Matrix<Rational> B = ones_vector<Rational>(PM.second.cols())|-T(PM.second);

      p.take("INEQUALITIES") << A/unit_matrix<Rational>(A.cols());
      Matrix<Rational> VP = p.give("VERTICES");

      Vector<Rational> za = unit_vector<Rational>(A.cols(),0);
      Vector<Rational> zb = unit_vector<Rational>(B.cols(),0);

      for ( Entire<Rows<Matrix<Rational> > >::const_iterator vit = entire(rows(VP)); !vit.at_end(); ++vit ) 
	if ( *vit != za ) {
	  Matrix<Rational> E = vector2row(0|A*(*vit));
	  E /= vector2row((*vit).slice(~scalar2set(0))*B);

	  perl::Object q("polytope::Polytope<Rational>");
	  q.take("INEQUALITIES") << B/unit_matrix<Rational>(B.cols());
	  q.take("EQUATIONS") << E;

	  if ( q.give("FEASIBLE") ) {
	    Matrix<Rational> VQ = q.give("VERTICES");
	    
	    std::pair<Vector<Rational>,Vector<Rational> > nash_equilibrium;
	    nash_equilibrium.second = (*vit).slice(~scalar2set(0))/(ones_vector<Rational>(vit->dim()-1)*(*vit).slice(~scalar2set(0)));

	    for ( Entire<Rows<Matrix<Rational> > >::const_iterator wit = entire(rows(VQ)); !wit.at_end(); ++wit ) 
	      if ( *wit != zb ) {
		nash_equilibrium.first = (*wit).slice(~scalar2set(0))/(ones_vector<Rational>(wit->dim()-1)*(*wit).slice(~scalar2set(0)));
		nash += nash_equilibrium;
	      }
	  }
      }
      return Array<std::pair<Vector<Rational>,Vector<Rational> > >(nash);
    }


    UserFunction4perl("# @category Game Theory", 
		      &nash_equilibria, "nash_equilibria(Pair<Matrix<Rational>,Matrix<Rational>>)");

  }}

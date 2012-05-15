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
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/Array.h"

namespace polymake { namespace bimatrix_games {

    // return the set of maximal entries of a vector
    Set<int> max_entries ( const Vector<Rational> & v ) {
      Set<int> s;
      Rational cur = v[0];  // the current maximal value
      s += 0;               // add its index to the set
      int i = 0;

      for ( Entire<Vector<Rational> >::const_iterator c = entire(v); !c.at_end(); ++c, ++i ) {
	if ( *c > cur ) { // reset maximum
	  cur = *c;
	  s.clear();
	  s += i;
	} else {          // add index
	  if ( *c == cur ) 
	    s += i;
	}
      }
      return s;
    }

    // find the support of a vector
    Set<int> support ( const Vector<Rational> & v ) {
      Set<int> s;
      int i = 0;
      for ( Entire<Vector<Rational> >::const_iterator c = entire(v); !c.at_end(); ++c, ++i ) 
	if ( *c != 0 ) s += i;
      return s;
    }

    // check whether a vector is either non-negative or non-positive
    bool sign_unique ( const Vector<Rational> & v ) {
      Vector<Rational>::const_iterator c = v.begin();
      while ( c != v.end() && *c == 0 ) ++c;
      const Rational first = *c;
      for ( Entire<Vector<Rational> >::const_iterator c = entire(v); !c.at_end(); ++c ) 
	if ( *c != 0 && sign(first) != sign(*c) ) return false;
      return true;
    }

    // add missing zeros to a vector, the index set of the entires in v is s, n is the final dimension
    Vector<Rational> extend ( const Vector<Rational> & v, const Set<int> & s, int n ) {
      Vector<Rational> vl(n);
      Vector<Rational>::const_iterator vit = v.begin();
      for ( Entire<Set<int> >::const_iterator sit = entire(s); !sit.at_end(); ++sit, ++vit ) 
	vl[*sit] = *vit;
      return vl;
    }


    // compute the nash equilibria of a non-generate game 
    // using the following algorithm:
    //
    Array<std::pair<Vector<Rational>,Vector<Rational> > > nash_equilibria_by_support_enumeration ( const std::pair<Matrix<Rational>,Matrix<Rational> > & PM ) {

      const Matrix<Rational>  A = PM.first;
      const Matrix<Rational>  B = PM.second;
     
      if ( A.rows() != B.rows() || A.cols() != B.cols() )   // not a bimatrix game
	throw std::runtime_error("dimension mismatch in input matrices");

      Set<std::pair<Vector<Rational>,Vector<Rational> > > nash;   // collects the nash equilibria
      int min = A.rows() < A.cols() ? A.rows() : A.cols();        // find the min of row and column number == max support of a mixed strategy
      
      
      for ( int k = 1; k <= min; ++k ) {
	Vector<Rational> ones = ones_vector<Rational>(k);
	for ( Entire<Subsets_of_k<const sequence& > >::const_iterator rit = entire(all_subsets_of_k(sequence(0,A.rows()),k)); !rit.at_end(); ++rit ) {
	  for ( Entire<Subsets_of_k<const sequence& > >::const_iterator cit = entire(all_subsets_of_k(sequence(0,A.cols()),k)); !cit.at_end(); ++cit ) {
	    Vector<Rational> x,y;
	    try {
	      x = lin_solve(T(B.minor(*rit,*cit)),ones);
	      y = lin_solve(A.minor(*rit,*cit),ones);
	    } catch (const std::exception& e) {
	      break;              // system might not have a solution
	                          // FIXME should distinguish between "no solution" and other errors 
	    }

	    if ( !sign_unique(x) || !sign_unique(y) ) break;   // not a nash equilibrium

	    Vector<Rational> xl = extend(x/(ones_vector<Rational>(x.dim())*x),*rit,A.rows());
	    Vector<Rational> yl = extend(y/(ones_vector<Rational>(y.dim())*y),*cit,A.cols());

	    if ( max_entries(A*yl) == support(xl)  &&  max_entries(T(B)*xl) == support(yl) ) 
	      nash += std::pair<Vector<Rational>,Vector<Rational> >(xl,yl);
	  }
	}
      }
      return Array<std::pair<Vector<Rational>,Vector<Rational> > >(nash);
    }


    UserFunction4perl("# @category Game Theory", 
		      &nash_equilibria_by_support_enumeration, "nash_equilibria_by_support_enumeration(Pair<Matrix<Rational>,Matrix<Rational>>)");

  }}

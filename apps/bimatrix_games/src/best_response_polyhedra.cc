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
#include "polymake/PowerSet.h"
#include "polymake/Set.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"

namespace polymake { namespace bimatrix_games {

    perl::Object first_best_response_polyhedron ( const Matrix<Rational> & A) {

      perl::Object q("polytope::Polytope<Rational>");

      Matrix<Rational> N(A);
      N = ones_vector<Rational>(N.rows())|-N;
      N = (zero_vector<Rational>(N.cols()-1)|unit_matrix<Rational>(N.cols()-1))/N;
      N = zero_vector<Rational>(N.rows())|N;

      q.take("INEQUALITIES") << N;
      q.take("EQUATIONS") << vector2row(-1|(0|ones_vector<Rational>(N.cols()-2)));

      return q;
    }


    perl::Object second_best_response_polyhedron ( const Matrix<Rational> & B) {

      perl::Object p("polytope::Polytope<Rational>");

      Matrix<Rational> M(T(B));
      M = ones_vector<Rational>(M.rows())|-M;
      M /= zero_vector<Rational>(M.cols()-1)|unit_matrix<Rational>(M.cols()-1);
      M = zero_vector<Rational>(M.rows())|M;

      p.take("INEQUALITIES") << M;
      p.take("EQUATIONS") << vector2row(-1|(0|ones_vector<Rational>(M.cols()-2)));

      return p;
    }

    perl::Object reduced_first_best_response_polyhedron ( const Matrix<Rational> & A ) {

      perl::Object q("polytope::Polytope<Rational>");

      Matrix<Rational> N(A);
      N = ones_vector<Rational>(N.rows())|-N;
      N = (zero_vector<Rational>(N.cols()-1)|unit_matrix<Rational>(N.cols()-1))/N;

      q.take("INEQUALITIES") << N;

      return q;
    }

    perl::Object reduced_second_best_response_polyhedron ( const Matrix<Rational> & B ) {

      perl::Object p("polytope::Polytope<Rational>");

      Matrix<Rational> M(T(B));
      M = ones_vector<Rational>(M.rows())|-M;
      M /= zero_vector<Rational>(M.cols()-1)|unit_matrix<Rational>(M.cols()-1);

      p.take("INEQUALITIES") << M;

      return p;
    }

    //    bool reduced_feasible ( const Matrix<Rational> & A, const Matrix<Rational> & B ) {
    bool reduced_feasible ( const std::pair<Matrix<Rational>,Matrix<Rational> > & PM ) {

      perl::Object p("polytope::Polytope<Rational>");
      perl::Object q("polytope::Polytope<Rational>");

      Matrix<Rational> M(T(PM.second));
      M = zero_vector<Rational>(M.rows())|-M;
      M /= zero_vector<Rational>(M.cols()-1)|unit_matrix<Rational>(M.cols()-1);

      p.take("INEQUALITIES") << M;
      int dp = p.give("CONE_DIM");

      Matrix<Rational> N(PM.first);
      N = zero_vector<Rational>(N.rows())|-N;
      N = (zero_vector<Rational>(N.cols()-1)|unit_matrix<Rational>(N.cols()-1))/N;

      q.take("INEQUALITIES") << N;
      int dq = q.give("CONE_DIM");

      return ( dp == 1 && dq == 1 );
    }


    UserFunction4perl("# @category Game Theory", 
		      &first_best_response_polyhedron, "first_best_response_polyhedron(Matrix)");

    UserFunction4perl("# @category Game Theory", 
		      &second_best_response_polyhedron, "second_best_response_polyhedron(Matrix)");

    UserFunction4perl("# @category Game Theory", 
		      &reduced_first_best_response_polyhedron, "reduced_first_best_response_polyhedron(Matrix)");

    UserFunction4perl("# @category Game Theory", 
		      &reduced_second_best_response_polyhedron, "reduced_second_best_response_polyhedron(Matrix)");

    UserFunction4perl("# @category Game Theory", 
		      &reduced_feasible, "reduced_feasible(Pair<Matrix<Rational>,Matrix<Rational>>)");


  }}

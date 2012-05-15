/* Copyright (c) 1997-2010
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
   $Project: polymake $$Id: wrap-connected_components.cc 224 2012-05-15 15:57:44Z paffenholz $
*/

namespace polymake { namespace bimatrix_games {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( pm::PowerSet<int, pm::operations::cmp> (pm::Array<std::pair<pm::Vector<pm::Rational>, pm::Vector<pm::Rational> >, void> const&) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Array< std::pair< Vector< Rational >, Vector< Rational > > > > >() );
   }
   FunctionWrapperInstance4perl( pm::PowerSet<int, pm::operations::cmp> (pm::Array<std::pair<pm::Vector<pm::Rational>, pm::Vector<pm::Rational> >, void> const&) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }

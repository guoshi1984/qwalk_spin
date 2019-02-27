/*
 
Copyright (C) 2007 Lucas K. Wagner

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
*/

#ifndef CCUBIC_SPLINE_H_INCLUDED
#define CCUBIC_SPLINE_H_INCLUDED

#include "CBasis_function.h"
#include "Spline_fitter.h"
#include "Real_Spherical_harmonics.h"
#include "Cubic_spline.h"
class Cubic_spline;


/*!
This class represents a 1-D cubic spline interpolation.  It is meant
to be used for the radial part of an atomic orbital, and will return
symmetry-adjusted values for the calculation methods.

*/


		     

class CCubic_spline: public CBasis_function
{ 
private:
  Cubic_spline * parent;

public:
  CCubic_spline()
  {}
  ~CCubic_spline()
  {}

  //-----------------------------------------------------

  int showinfo(string & indent, ostream & os);
  
  int writeinput(string &, ostream &);
  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos
    //!< The current position in the words(important if one basis section makes several functions); will be incremented as the Basis_function reads the words.
  );

  int nfunc();
  doublevar cutoff(int);
  
  void raw_input(ifstream & input);
  
  string label();

  void calcVal(const Array1 <doublevar> & r,
               Array1 <dcomplex> & symvals,
               const int startfill=0);

  /*!
     All returned values are of the form \f$f(r),
     \frac{1}{r}\frac{df(r)}{dr}, \frac{d^2f(r)}{dr^2} \f$.
  */
  void calcLap(
    const Array1 <doublevar> & r,
    //!< in form r, r^2, x, y, z
    Array2 <dcomplex> & symvals,
    //!< The values of the spline propogated through symmetry.  For example, a p state would be a 3x5 matrix, and an s state a 1x5.
    const int startfill=0
  );

 };

#endif // CCUBIC_SPLINE_H_INCLUDED
//--------------------------------------------------------------------------

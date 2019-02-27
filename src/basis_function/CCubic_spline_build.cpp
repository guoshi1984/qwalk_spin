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

#include "Qmc_std.h"
#include "qmc_io.h"
#include <iomanip>
#include "CCubic_spline.h"

//-------------------------------------------------------------------------
int CCubic_spline::read(
  vector <string> & words,
  unsigned int & pos
)
{
  parent=new Cubic_spline;
  return parent->read(words, pos);
}
//-------------------------------------------------------------------------
void CCubic_spline::raw_input(ifstream & input)
{
  parent->raw_input(input);
}
//-------------------------------------------------------------
int CCubic_spline::nfunc()
{
  return parent->nfunc();
}
//----------------------------------------------------------------------

doublevar CCubic_spline::cutoff(int n){
  return parent->cutoff(n);
}
//----------------------------------------------------------------------
string CCubic_spline::label()
  {
    return parent->label();
  }

//----------------------------------------------------------------------
int CCubic_spline::showinfo(string & indent, ostream & os)
{
  os << indent << "Cubic spline for " << parent->atomname << endl;
  if(parent->zero_derivative)
    os << indent << "Zero derivative enforced at the orgin\n";

  if(parent->recursive){
    os << indent << "Using recursive evaluation of complex spherical harmonics (addapted from Jeognim Kim's QMCPACK) upto Lmax "<<parent->Lmax<<endl;
  }

  if(parent->customspacing!=0.02)
    os << indent << "Using custom spacing of "<<parent->customspacing<<endl;
  os << indent << parent->nsplines << "  radial functions\n";
  os << indent << setw(10) << "function" << setw(10) << "symmetry"
     << setw(10) << "cutoff" << endl;
  
  

  
  int totfunc=0;
  for(int i=0; i< parent->nsplines; i++) {
    os << indent << setw(10) << i << setw(10) << parent->symmetry_lookup(parent->symmetry(i))
       << setw(10) << parent->rcut(totfunc) << endl;
    
    totfunc+=parent->nfuncspline(i);
    
  }

  return 1;
}

//----------------------------------------------------------------------

int CCubic_spline::writeinput(string & indent, ostream & os)
{
  return parent->writeinput(indent, os);
}
//------------------------------------------------------------------------

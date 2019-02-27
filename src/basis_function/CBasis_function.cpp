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


#include "CBasis_function.h"
#include "qmc_io.h"
#include "CPlanewave_function.h"
#include "Blochwave_function.h"
#include "CCubic_spline.h"
#include "SO_CPlanewave_function.h"
#include <iostream>

int allocate(vector <string> & basistext, CBasis_function * & bptr)
{ 
  assert( bptr == NULL );
  unsigned int pos=0;
  string type;
  readnext(basistext, pos, type);
  if(caseless_eq(type,"CPLANEWAVE"))
    bptr=new CPlanewave_function;
  else if(caseless_eq(type,"BLOCHWAVE"))
    bptr=new Blochwave_function;
  else if(caseless_eq(type ,"AOSPLINE"))
    bptr= new CCubic_spline;
  else if(caseless_eq(type,"SO_CPLANEWAVE"))
  {  bptr= new SO_CPlanewave_function;  }
  else
    error("Didn't understand the basis section type ", type,
          " in basis section ", basistext[0]);
   
  bptr->read(basistext, pos);
  return 1;

}

int deallocate(CBasis_function * & bptr)
{
  if(bptr == NULL)
    return 0;

  delete bptr;
  bptr=NULL;
  return 1;
}
//--------------------------------------------------------------------------

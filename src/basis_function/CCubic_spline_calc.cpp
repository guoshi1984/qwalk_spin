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
#include <iomanip>
#include "CCubic_spline.h"
#include "Cubic_spline.h"


void CCubic_spline::calcVal(const Array1 <doublevar> & r,
                           Array1 <dcomplex> & symvals,
                           const int startfill)
{
  const dcomplex I(0.0,1.0);
  const doublevar sqrt2 = sqrt(2.0);
  if(r(0) >= parent->threshold)
    {
    int end=startfill+parent->nfunctions;
    for(int i=startfill; i< end; i++)
    {
      symvals(i)=0;
    }
  }
  else {
    if(!parent->recursive){ //original approach
      error("only the recursive evaluation is implemented for CCubic_spline");
    }
    else { //recursive evaluation
      parent->realYlm->evaluate(r);
      doublevar temp;
      int totf=startfill;
      for(int i=0; i<parent->nsplines; i++)  {
	int interval=parent->splines(i).getInterval(r(0));
	temp=parent->splines(i).getVal(r(0), interval);
	int l=parent->symmetry_lvalue(parent->symmetry(i));
	doublevar factor; //(-1)^m \sqrt{2}
	for(int m=l;m>=-l;m--){
	  factor=1.0;
	  if(m==0)
	    symvals(totf++)=temp*parent->realYlm->getYlm(l, m); //only real part 
	  else if (m>0){
	    if(parent->addsign)
	      factor=pow(-1.0e0,m)/sqrt2;
	    else
	      factor=1.0/sqrt2;

	    symvals(totf++)=factor*temp*(parent->realYlm->getYlm(l, m) + I*parent->realYlm->getYlm(l, -m)); 
	  }
	  else{ //m<0
	    if(parent->addsign)
//             Duplicate dividing by factor sqrt(2.0) corrected 		     
//	       factor=pow(-1.0e0,-m)/sqrt2;
	       factor=pow(-1.0e0,-m);
	     else
//	       factor=1.0/sqrt2;
	       factor=1.0;	     
//	    symvals(totf++)=factor*conj(symvals(totf-2*abs(m)-1));
	    symvals(totf)=factor*conj(symvals(totf-2*abs(m)));
	    totf++;
	  }
	}
      }//end of loop over splines
    }//end of recursive
  }
}

void CCubic_spline::calcLap(
  const Array1 <doublevar> & r,
  //!< in form r, r^2, x, y, z
  Array2 <dcomplex> & symvals,
  //!< The values of the spline propogated through symmetry.  For example, a p state would be a 3x3 matrix, and an s state a 1x3.
  const int startfill){

  const dcomplex I(0.0,1.0);
  const doublevar sqrt2 = sqrt(2.0);
  const doublevar two=2.0;
  
  assert(r.GetDim(0) >= 5);
  //cout << "spline interval " << interval << "   r   " << r(0) << endl;
  if(r(0) >= parent->threshold)
  {
    int end=startfill+parent->nfunctions;
    for(int i=startfill; i< end; i++)
    {
      for(int d=0; d< 5; d++)
      {
        symvals(i,d)=0;
      }
    }
  }
  else
  {
    if(!parent->recursive){ //original approach
      error("only the recursive evaluation is implemented for CCubic_spline");
    }
    else { //recursive approach 
      register doublevar func, fdir, f2dir;
      register doublevar x, y, z;
      register dcomplex h, v4, dfdx, dfdy, dfdz;
      assert(r(0) < parent->threshold);
      assert(symvals.GetDim(0) >= parent->nfunctions);

       parent->realYlm->evaluateAll(r);
       Array1 <doublevar> gYlm(3);

       x=r(2);
       y=r(3);
       z=r(4);

       int totf=startfill;
       doublevar factor;
       for(int i=0; i<parent->nsplines; i++)  {
	 int interval=parent->splines(i).getInterval(r(0));
	 parent->splines(i).getDers(r(0), interval,func, fdir, f2dir);
	 int l=parent->symmetry_lvalue(parent->symmetry(i));
	 for(int m=l;m>=-l;m--){
	   if(m==0){ //only real part 
	      v4=parent->realYlm->getYlm(l, m);
	      parent->realYlm->getGradYlm(l,m,gYlm);
	      dfdx=gYlm(0);
	      dfdy=gYlm(1);
	      dfdz=gYlm(2);
	      h=fdir*v4;
	      symvals(totf,0)=v4*func;
	      symvals(totf,1)=func*dfdx+h*x;
	      symvals(totf,2)=func*dfdy+h*y;
	      symvals(totf,3)=func*dfdz+h*z;
	      symvals(totf,4)=two*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	   }
	   else if (m>0){
	     if(parent->addsign)
	      factor=pow(-1.0e0,m)/sqrt2;
	     else
	      factor=1.0/sqrt2;

	     v4=factor*(parent->realYlm->getYlm(l, m)+I*parent->realYlm->getYlm(l, -m));
	     parent->realYlm->getGradYlm(l,m,gYlm);
	     dfdx=gYlm(0);
	     dfdy=gYlm(1);
	     dfdz=gYlm(2);

	     parent->realYlm->getGradYlm(l,-m,gYlm);
	     dfdx+=I*gYlm(0);
	     dfdy+=I*gYlm(1);
	     dfdz+=I*gYlm(2);

	     dfdx*=factor;
	     dfdy*=factor;
	     dfdz*=factor;

	     h=fdir*v4;
	     symvals(totf,0)=v4*func;
	     symvals(totf,1)=func*dfdx+h*x;
	     symvals(totf,2)=func*dfdy+h*y;
	     symvals(totf,3)=func*dfdz+h*z;
	     symvals(totf,4)=two*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	   }
	   else{ //m<0
	     if(parent->addsign)
//             Duplicate dividing by factor sqrt(2.0) corrected 		     
//	       factor=pow(-1.0e0,-m)/sqrt2;
	       factor=pow(-1.0e0,-m);
	     else
//	       factor=1.0/sqrt2;
	       factor=1.0;
	     symvals(totf,0)=factor*conj(symvals(totf-2*abs(m), 0));
	     symvals(totf,1)=factor*conj(symvals(totf-2*abs(m), 1));
	     symvals(totf,2)=factor*conj(symvals(totf-2*abs(m), 2));
	     symvals(totf,3)=factor*conj(symvals(totf-2*abs(m), 3));
	     symvals(totf,4)=factor*conj(symvals(totf-2*abs(m), 4));
	   }
	   totf++;
	 }
       }//end of loop over splines
    }//end of recursive
  }
}



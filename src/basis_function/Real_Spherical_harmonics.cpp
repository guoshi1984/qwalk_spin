//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim, John Shumway and J. Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
//Addapted from QMCPACK's SphericalTensor class with agreement from Jeongnim Kim. M.B.

#include "Real_Spherical_harmonics.h"

void RealSphericalTensor::init(const int l_max, bool addsign){ 
  Lmax=l_max;
  int ntot = (Lmax+1)*(Lmax+1);
  Ylm.Resize(ntot);
  gradYlm.Resize(ntot);
  HessYlm.Resize(ntot);
  for(int i=0;i<ntot;i++){
    gradYlm(i).Resize(3);
    gradYlm(i) = 0.0;
    HessYlm(i).Resize(3,3);
    HessYlm(i) = 0.0;
  }
  

  NormFactor.Resize(ntot);
  NormFactor=1.0;
  const doublevar sqrt2 = sqrt(2.0);
  if(addsign) {
    for (int l=0; l<=Lmax; l++) {
      NormFactor(index(l,0))=1.0;
      for (int m=1; m<=l; m++) {
        NormFactor(index(l,m))=pow(-1.0e0,m)*sqrt2;
        NormFactor(index(l,-m))=pow(-1.0e0,-m)*sqrt2;
     }
    }
  } else {
    for (int l=0; l<=Lmax; l++) {
      for (int m=1; m<=l; m++) {
	NormFactor(index(l,m))=sqrt2;
	NormFactor(index(l,-m))=sqrt2;
     }
    }
  }

  FactorL.Resize(Lmax+1);
  const doublevar omega = 1.0/sqrt(16.0*atan(1.0));
  for(int l=1; l<=Lmax; l++) FactorL(l) = sqrt(static_cast<doublevar>(2*l+1))*omega;

  Factor2L.Resize(Lmax+1);
  for(int l=1; l<=Lmax; l++) Factor2L(l) = static_cast<doublevar>(2*l+1)/static_cast<doublevar>(2*l-1);

  FactorLM.Resize(ntot);
  for(int l=1; l<=Lmax; l++) 
    for(int m=1; m<=l; m++) {
      doublevar fac2 = 1.0/sqrt(static_cast<doublevar>((l+m)*(l+1-m)));
      FactorLM(index(l,m))=fac2;
      FactorLM(index(l,-m))=fac2;
    }
}

void RealSphericalTensor::evaluate(const Array1 <doublevar> & p) {
  value_type x=p(2), y=p(3), z=p(4);
  const value_type pi = 4.0*atan(1.0);
  const value_type pi4 = 4.0*pi;
  const value_type omega = 1.0/sqrt(pi4);
  const value_type sqrt2 = sqrt(2.0);

  /*  Calculate r, cos(theta), sin(theta), cos(phi), sin(phi) from input 
      coordinates. Check here the coordinate singularity at cos(theta) = +-1. 
      This also takes care of r=0 case. */

  value_type cphi,sphi,ctheta;
  value_type r2xy=x*x+y*y;
  value_type r=sqrt(r2xy+z*z);  
  if (r2xy<numeric_limits<T>::epsilon()) {
     cphi = 0.0;
     sphi = 1.0;
     ctheta = (z<0)?-1.0:1.0;
  } else {
     ctheta = z/r;
     value_type rxyi = 1.0/sqrt(r2xy);
     cphi = x*rxyi;
     sphi = y*rxyi;
  }

  value_type stheta = sqrt(1.0-ctheta*ctheta);

  /* Now to calculate the associated legendre functions P_lm from the
     recursion relation from l=0 to Lmax. Conventions of J.D. Jackson, 
     Classical Electrodynamics are used. */
  
  Ylm(0) = 1.0;

  // calculate P_ll and P_l,l-1

  value_type fac = 1.0; 
  int j = -1;
  for (int l=1; l<=Lmax; l++) {
    j += 2;
    fac *= -j*stheta;
    int ll=index(l,l);
    int l1=index(l,l-1);
    int l2=index(l-1,l-1);
    Ylm(ll) = fac;
    Ylm(l1) = j*ctheta*Ylm(l2);
  }

  // Use recurence to get other plm's //
  for (int m=0; m<Lmax-1; m++) {
    int j = 2*m+1;
    for (int l=m+2; l<=Lmax; l++) {
      j += 2;
      int lm=index(l,m);
      int l1=index(l-1,m);
      int l2=index(l-2,m);
      Ylm(lm) = (ctheta*j*Ylm(l1)-(l+m-1)*Ylm(l2))/(l-m);
    }
  }
  
  // Now to calculate r^l Y_lm. //
  value_type sphim,cphim,fac2,temp;
  Ylm(0) = omega; //1.0/sqrt(pi4);          
  value_type rpow = 1.0;
  for (int l=1; l<=Lmax; l++) {
    rpow *= r;
    //fac = rpow*sqrt(static_cast<T>(2*l+1))*omega;//rpow*sqrt((2*l+1)/pi4);  
    //FactorL[l] = sqrt(2*l+1)/sqrt(4*pi)
    fac = rpow*FactorL(l);
    int l0=index(l,0);
    Ylm(l0) *= fac;
    cphim = 1.0;
    sphim = 0.0;
    for (int m=1; m<=l; m++) {
      temp = cphim*cphi-sphim*sphi;
      sphim = sphim*cphi+cphim*sphi;
      cphim = temp;
      int lm = index(l,m);    
      fac *= FactorLM(lm);
      temp = fac*Ylm(lm);
      Ylm(lm) = temp*cphim;
      lm = index(l,-m);
      Ylm(lm) = temp*sphim;
    }
  }
  for (int i=0; i<Ylm.GetSize(); i++) Ylm(i)*= NormFactor[i];
}

void RealSphericalTensor::evaluateAll(const Array1 <doublevar> & p) {
  value_type x=p[2], y=p[3], z=p[4];
  const value_type pi = 4.0*atan(1.0);
  const value_type pi4 = 4.0*pi;
  const value_type omega = 1.0/sqrt(pi4);
  const value_type sqrt2 = sqrt(2.0);

  //cout << " <numeric_limits<T>::epsilon() = "<< numeric_limits<T>::epsilon() <<endl;

  /*  Calculate r, cos(theta), sin(theta), cos(phi), sin(phi) from input 
      coordinates. Check here the coordinate singularity at cos(theta) = +-1. 
      This also takes care of r=0 case. */

  value_type cphi,sphi,ctheta;
  value_type r2xy=x*x+y*y;
  value_type r=sqrt(r2xy+z*z);  
  if (r2xy<numeric_limits<T>::epsilon()) {
     cphi = 0.0;
     sphi = 1.0;
     ctheta = (z<0)?-1.0:1.0;
  } else {
     ctheta = z/r;
     value_type rxyi = 1.0/sqrt(r2xy);
     cphi = x*rxyi;
     sphi = y*rxyi;
  }

  value_type stheta = sqrt(1.0-ctheta*ctheta);

  /* Now to calculate the associated legendre functions P_lm from the
     recursion relation from l=0 to Lmax. Conventions of J.D. Jackson, 
     Classical Electrodynamics are used. */
  Ylm(0) = 1.0;

  // calculate P_ll and P_l,l-1

  value_type fac = 1.0; 
  int j = -1;
  for (int l=1; l<=Lmax; l++) {
    j += 2;
    fac *= -j*stheta;
    int ll=index(l,l);
    int l1=index(l,l-1);
    int l2=index(l-1,l-1);
    Ylm(ll) = fac;
    Ylm(l1) = j*ctheta*Ylm(l2);
  }

  // Use recurence to get other plm's //
  
  for (int m=0; m<Lmax-1; m++) {
    int j = 2*m+1;
    for (int l=m+2; l<=Lmax; l++) {
      j += 2;
      int lm=index(l,m);
      int l1=index(l-1,m);
      int l2=index(l-2,m);
      Ylm(lm) = (ctheta*j*Ylm(l1)-(l+m-1)*Ylm(l2))/(l-m);
    }
  }
  
  // Now to calculate r^l Y_lm. //
  
  value_type sphim,cphim,fac2,temp;
  Ylm[0] = omega; //1.0/sqrt(pi4);          
  value_type rpow = 1.0;
  for (int l=1; l<=Lmax; l++) {
    rpow *= r;
    //fac = rpow*sqrt(static_cast<T>(2*l+1))*omega;//rpow*sqrt((2*l+1)/pi4);  
    //FactorL[l] = sqrt(2*l+1)/sqrt(4*pi)
    fac = rpow*FactorL(l);
    int l0=index(l,0);
    Ylm(l0) *= fac;
    cphim = 1.0;
    sphim = 0.0;
    for (int m=1; m<=l; m++) {
      //fac2 = (l+m)*(l+1-m);
      //fac = fac/sqrt(fac2);
      temp = cphim*cphi-sphim*sphi;
      sphim = sphim*cphi+cphim*sphi;
      cphim = temp;
      int lm = index(l,m);    
      //fac2,fac use precalculated FactorLM
      fac *= FactorLM(lm);
      temp = fac*Ylm(lm);
      Ylm(lm) = temp*cphim;
      lm = index(l,-m);
      Ylm(lm) = temp*sphim;
    }
  }

  // Calculating Gradient now//

  for (int l=1; l<=Lmax; l++) {
    //value_type fac = ((value_type) (2*l+1))/(2*l-1);
    value_type fac = Factor2L(l);
    for (int m=-l; m<=l; m++) {
      int lm = index(l-1,0);  
      value_type gx,gy,gz,dpr,dpi,dmr,dmi;
      int ma = abs(m);
      value_type cp = sqrt(fac*(l-ma-1)*(l-ma));
      value_type cm = sqrt(fac*(l+ma-1)*(l+ma));
      value_type c0 = sqrt(fac*(l-ma)*(l+ma));
      gz = (l > ma) ? c0*Ylm(lm+m):0.0;
      if (l > ma+1) {
        dpr = cp*Ylm(lm+ma+1);
        dpi = cp*Ylm(lm-ma-1);
      } else {
        dpr = 0.0;
        dpi = 0.0;
      }
      if (l > 1) {
        switch (ma) {
          
        case 0:
          dmr = -cm*Ylm(lm+1);
          dmi = cm*Ylm(lm-1);
          break;
        case 1:
          dmr = cm*Ylm(lm);
          dmi = 0.0;
          break;
        default:
          dmr = cm*Ylm(lm+ma-1);
          dmi = cm*Ylm(lm-ma+1);
        }
      } else {
        dmr = cm*Ylm(lm);
        dmi = 0.0;
        //dmr = (l==1) ? cm*Ylm[lm]:0.0;
        //dmi = 0.0;
      }
      if (m < 0) {
        gx = 0.5*(dpi-dmi);
        gy = -0.5*(dpr+dmr);
      } else {
        gx = 0.5*(dpr-dmr);
        gy = 0.5*(dpi+dmi);
      }
      lm = index(l,m);
      Array1 <doublevar> gtmp(3);
      
      if(ma){
	gtmp(0)=NormFactor(lm)*gx;
	gtmp(1)=NormFactor(lm)*gy;
	gtmp(2)=NormFactor(lm)*gz;
        gradYlm(lm)=gtmp; 
      }
      else{
	gtmp(0)=gx;
	gtmp(1)=gy;
	gtmp(2)=gz;
        gradYlm(lm) = gtmp; 
      }
    }
  }
  for (int i=0; i<Ylm.GetSize(); i++) Ylm(i)*= NormFactor(i);
 //for (int i=0; i<Ylm.size(); i++) gradYlm[i]*= NormFactor[i];
}


void RealSphericalTensor::evaluateincludingHessian(const Array1 <doublevar> & p) {
  error("Recursive evaluation of Hessian in RealSphericalTensor not implemented yet so backflow will not work!, M.B.");
}

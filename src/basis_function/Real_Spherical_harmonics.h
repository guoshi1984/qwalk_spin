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


#ifndef REAL_SPHERICAL_HARMONICS_H_INCLUDED
#define REAL_SPHERICAL_HARMONICS_H_INCLUDED
#include "CBasis_function.h"
#include <limits>

class RealSphericalTensor {
public : 
  typedef doublevar value_type;
  typedef value_type T;
  
  /** constructor
   * @param l_max maximum angular momentum
   * @param addsign flag to determine what convention to use
   *
   * Evaluate all the constants and prefactors.
   * The spherical harmonics is defined as
   * \f[ Y_l^m (\theta,\phi) = \sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}} P_l^m(\cos\theta)e^{im\phi}\f]
   * Note that the data member Ylm is a misnomer and should not be confused with "spherical harmonics" 
   * \f$Y_l^m\f$.
   - When addsign == true, e.g., Gaussian packages
   \f{eqnarray*}
     S_l^m &=& (-1)^m \sqrt{2}\Re(Y_l^{|m|}), \;\;\;m > 0 \\
     &=& Y_l^0, \;\;\;m = 0 \\
     &=& (-1)^m \sqrt{2}\Im(Y_l^{|m|}),\;\;\;m < 0
   \f}
   - When addsign == false, e.g., SIESTA package,
   \f{eqnarray*}
     S_l^m &=& \sqrt{2}\Re(Y_l^{|m|}), \;\;\;m > 0 \\
     &=& Y_l^0, \;\;\;m = 0 \\
     &=&\sqrt{2}\Im(Y_l^{|m|}),\;\;\;m < 0
   \f}
  */

  RealSphericalTensor(){}
  ~RealSphericalTensor(){}

  void init(const int l_max, bool addsign);

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluate(const Array1 <doublevar> & r);     

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateAll(const Array1 <doublevar> & r);     

  void evaluateincludingHessian(const Array1 <doublevar> & r);     

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateTest(const Array1 <doublevar> & r);     

  ///returns the index/locator for (\f$l,m\f$) combo, \f$ l(l+1)+m \f$
  inline int index(int l, int m) const {return (l*(l+1))+m;}

  ///returns the value of \f$ r^l S_l^m \f$ given \f$ (l,m) \f$
  doublevar getYlm(int l, int m) const
  {return Ylm(index(l,m));}

  ///returns the gradient of \f$ r^l S_l^m \f$ given \f$ (l,m) \f$
  void getGradYlm(int l, int m, Array1 < doublevar> & gYlm)
  { gYlm=gradYlm(index(l,m)); }

  ///returns the value of \f$ r^l S_l^m \f$ given \f$ (l,m) \f$
  doublevar getYlm(int lm) const {return Ylm(lm);}

  ///returns the gradient of \f$ r^l S_l^m \f$ given \f$ (l,m) \f$
  void getGradYlm(int lm, Array1 < doublevar> & gYlm) 
  { gYlm=gradYlm(lm); }
  
  void getHessYlm(int lm, Array2 < doublevar> & hYlm) 
  { hYlm=HessYlm(lm); }

  void getHessYlm(int l, int m, Array2 < doublevar> & hYlm)
  { hYlm=HessYlm(index(l,m)); }


  inline int size() const { return Ylm.GetSize();}

  inline int lmax() const { return Lmax;}


  ///maximum angular momentum for the center
  int Lmax;

private:
  ///values  Ylm\f$=r^l S_l^m(x,y,z)\f$
  Array1 <doublevar> Ylm;
  /// Normalization factors 
  Array1 <doublevar> NormFactor;
  ///pre-evaluated factor \f$1/\sqrt{(l+m)\times(l+1-m)}\f$
  Array1 <doublevar> FactorLM;
  ///pre-evaluated factor \f$\sqrt{(2l+1)/(4\pi)}\f$
  Array1 <doublevar> FactorL;
  ///pre-evaluated factor \f$(2l+1)/(2l-1)\f$
  Array1 <doublevar> Factor2L;
  ///gradients gradYlm\f$=\nabla r^l S_l^m(x,y,z)\f$
  Array1 < Array1 < doublevar> > gradYlm;
  Array1 < Array2 < doublevar> > HessYlm;

};

#endif

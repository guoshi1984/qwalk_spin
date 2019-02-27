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

#ifndef SO_CSLAT_WF_H_INCLUDED

#define SO_CSLAT_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
class Wavefunction_data;
class Slat_wf_data;
class System;


class SO_Cslat_wf_storage : public Wavefunction_storage
{
public:
  virtual ~SO_Cslat_wf_storage()
  {}
private:
  friend class SO_Cslat_wf;

  //dimensions are [value gradient lap, MO]
  //increase the dimension to save the spinor
  Array2 <Array1 <dcomplex> >  moVal_temp;
  //should change the dimension
  Array3 <Array2 <dcomplex> >  inverse_temp;
  Array3 <dcomplex> detVal_temp;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
class SO_Cslat_wf : public  Wavefunction
{

public:

  SO_Cslat_wf()
  {}


  virtual int nfunc() {
    return nfunc_;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);
  virtual void updateSpinOrbit(Wavefunction_data *, Sample_point *);
  virtual void updateSpinProj(Wavefunction_data *, Sample_point* ); 
  virtual void getSpin(Wavefunction_data*, int, Array2<dcomplex>& value ); 
  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);
  virtual void getSpinOrbit(Wavefunction_data*, int, Array1<dcomplex>& );
  virtual void getDensity(Wavefunction_data *,int, Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return &);

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  void generateStorage(Wavefunction_storage * & wfstore);

  
  void init(Wavefunction_data *);

//  void getSpinProjVal(Wavefunction_data *, int, Wf_return &);
  
  
  //--
private:

  void save_for_static();

  void calcVal(Sample_point *);
  void updateVal(Sample_point *, int);
  void calcLap(Sample_point *);
  void updateLap(Sample_point *, int);

  void calcSpinOrbit(Sample_point *);
  void updateSpinOrbit(Sample_point *, int);

  void calcSpinProj(Sample_point* sample);
  void updateSpinProj(Sample_point* sample,int);

  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  Array1 <int> electronIsStaleSO;
  Array1 <int> electronIsStaleSpin;

  int updateEverythingVal;
  int updateEverythingLap;
  int updateEverythingSO;
  int updateEverythingSpin;

  int sampleAttached;
  int dataAttached;
  Slat_wf_data * parent;

  //Saved variables for electron updates, 5 components, elect, ith mo
  Array3 < Array1 <dcomplex> >  moVal;

  //change dimensions to 3 mo_index, 5 components, spinor
  Array3 <dcomplex> updatedMoVal;

  //Array2 <dcomplex> updatedSpinProjMoVal;

  Array3 < Array2 <dcomplex> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array3 <dcomplex> detVal;

//  Array3 <dcomplex> SpinProjMoVal; //proj index, electron index, mo index
  //Variables for a static(electrons not moving) calculation

  Array2 <dcomplex> updatedMoSO;
  Array3 <dcomplex> MoSOVal;

  Array2 <dcomplex> updatedSpinProjMoVal;
  Array3 <dcomplex> SpinProjMoVal;
  int staticSample;
  Array3 <dcomplex> saved_laplacian;
  //!<Saved laplacian for a static calculation (electron, function, [val grad lap])


  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron

};

#endif //SO_CSLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------

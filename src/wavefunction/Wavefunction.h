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

#ifndef WAVEFUNCTION_H_INCLUDED
#define WAVEFUNCTION_H_INCLUDED

#include "qmc_io.h"
#include "Qmc_std.h"
#include "Array.h"
class System;
class Sample_point;
class Wavefunction_data;
#include <vector>
#include <string>
#include <cstdio>
#include <complex>

/*!
Let a wavefunction know what changes have happened.
*/
enum change_type {sample_attach,//!< A sample point is reporting to this wf
                  data_attach,  //!< A wf_data is reporting
                  electron_move,//!< One electron changed position
                  all_electrons_move,//!< All electrons changed position
                  wf_parm_change,//!< One wavefunction parameter changed
                  all_wf_parms_change,//!< all wavefunction parameters changed
                  ion_move,//!< One ion(atom) moved
                  all_ions_move,//!< all atoms moved
                  sample_static,//!< We promise that we won't move an electron
                  sample_dynamic//!< revoke the promise of sample_static
                  };


/*!
Stores the state of the wavefunction.
*/
class Wavefunction_storage
{
public:
  virtual ~Wavefunction_storage()
  {}

private:
};

typedef complex <doublevar> dcomplex;
//typedef Array2 <complex <doublevar> > Wf_return;


struct Wf_return {
  Wf_return() {is_complex=0;}
  
  void write(string & indent, ostream & os);
  void read(istream & is);
  void mpiSend(int node);
  void mpiRecieve(int node);
  
  Wf_return(int nfunc, int nst) {
    is_complex=0;
    amp.Resize(nfunc, nst);
    phase.Resize(nfunc, nst);
    cvals.Resize(nfunc, nst);
  }
  void Resize(int nfunc, int nst) {
    amp.Resize(nfunc, nst);
    phase.Resize(nfunc, nst);
    cvals.Resize(nfunc, nst);
  }

  int sign(int w) {
    if(is_complex) return 1;

    if(fabs(phase(w,0)) < 1e-6) return 1;
    else return -1;
  }


  /*!
    used for complex functions
    vals= [ln|psi|, grad psi/psi, grad^2 psi/psi
    p=phase
   */
  void setVals(Array2 <dcomplex> & vals, Array1 <doublevar> & p);

  /*!
    used for real functions
    vals= [ln|psi|, grad psi/psi, grad^2 psi/psi
    sign= sign of wave function
  */
  void setVals(Array2 <doublevar> & vals, Array1 <doublevar> &  sign);

  int is_complex;
  Array2 <doublevar> amp;//!< ln( |psi| ), grad ln( |psi| ), grad^2 |psi|/|psi|
  Array2 <doublevar> phase; //!< phase and derivatives
  Array2 <dcomplex> cvals; //!< (null), grad ln(|psi|), grad^2 psi/psi  for debuggin purposes.
  Array2 <dcomplex> spinproj;
};



struct Parm_deriv_return { 
  int need_hessian;
  int nparms_start,nparms_end;
  int is_complex;
  Array1 <doublevar> gradient;
  Array2 <doublevar> hessian;
  Array1 <dcomplex> cgradient;
  Array2 <dcomplex> chessian;
  Parm_deriv_return() { 
    need_hessian=0;
    nparms_start=nparms_end=0;
    is_complex=0;
  }
};

/*!  extend the Hessian matrix assuming that the variables are independent
*/
void extend_parm_deriv(Parm_deriv_return & ret1, const Parm_deriv_return & ret2);


/*!
A base class for wavefunctions..

The basic functions are:

notify - notifies the wavefunction that something has
changed.  Usually used only by Sample_point and Wavefunction_data,
since they are the observed quantities.

updateprop - Brings the state of the wavefunction up to date
with the things that it's observing.

getprop - Returns the last property values.

where prop can be any of Val,  Lap, or ForceBias

Val gives you only the value

ForceBias is for vmc type methods, where you need importance sampling, but
also want a fast evaluation, so ForceBias is a relatively fast approximation
to the value and gradients.  It is not
really used any more, though, and defaults to the Lap methods unless 
overridden by the child class.

Lap gives you all the quantities, exactly.
*/
class Wavefunction
{
private:
protected:

public:
  virtual ~Wavefunction()
  {}


  Wavefunction()
  {}


  /*!
    \brief
    Notify the object of a change.  @see change_type
   */
  virtual void notify(change_type , int )=0;

  /*!
    number of functions that this object represents(ie, the length of the
    first index in the getX returns)
   */
  virtual int nfunc() {
    return 1;
  }

  /*!
    \brief
    Initialize with a wavefunction_data(should it be private, since
    wf_data is a friend, usually?)
  */
  virtual void init(Wavefunction_data * )
  {error("This Wavefunction object doesn't have init");}


  virtual void updateVal(Wavefunction_data *, Sample_point *)=0;
  virtual void updateLap(Wavefunction_data *, Sample_point *)=0;
 
  virtual void updateForceBias(Wavefunction_data *, Sample_point* );

  virtual void updateSpinProj(Wavefunction_data *, Sample_point* );
  virtual void updateSpinOrbit(Wavefunction_data *, Sample_point *);
 
  virtual void getVal(Wavefunction_data *, int, Wf_return &)=0;

  virtual void getLap(Wavefunction_data *, int, Wf_return &)=0;
  
  virtual void getForceBias(Wavefunction_data *, int, Wf_return &);

  virtual void getSpinOrbit(Wavefunction_data *, int, Array1<dcomplex> &);
                                                 // array for every wf
  virtual void getSpin(Wavefunction_data*, int, Array2<dcomplex>& );
  /*!
    \brief
    Calculate the derivatives with respect to the parameters of
    the wave function.  Does not change the state of the wave function,
    other than to update values

    Will resize gradient and hessian in Parm_deriv_return to 
    the correct size.
  */
  virtual int getParmDeriv(Wavefunction_data *, Sample_point *,
			    Parm_deriv_return & ) {
    return 0;
  }

  /*!
    \brief
    return symetric part of WF value (needed to get only Jastrow part in SH-DMC)
   */
  virtual void getSymmetricVal(Wavefunction_data *, int, Wf_return &)=0;


  /*!
    \brief
    get an approximation to the one-particle density(use updateVal to update)
   */
  virtual void getDensity(Wavefunction_data *,int,  Array2 <doublevar> &)=0;

  /*!
    \brief
    generate the correct Wavefunction_storage object to store an electron
    update
   */
  virtual void generateStorage(Wavefunction_storage * & wfstore)=0;


  /*!
    \brief
    saves an electron update to the storage object
   */
  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *)=0;
  
  /*!
    \brief
    saves an two electron update to the storage object
    */
  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *)
  {error("This Wavefunction object doesn't have two electron storage");}
   
  /*!
    \brief
    restores an electron update from the storage object(previously stored      with saveUpdate).

    Note that it only has a memory of one electron.  If you store,
    move two or more electrons, and then restore, the behavior is undefined.
   */
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *)=0;
  
  /*!
   \brief
   Restores  situation after two electron update
   */
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *)
  {error("This Wavefunction object doesn't have two electron storage");}


  /*!
    \brief
    Stores wavefunction parameter independent values to the given array.

    Find out the size of the array from Wavefunction_data::valSize().
    There are no assertions over what is in the array, or what order,
    so it should be treated as a black box.
   */
  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & )=0;
  /*!
    \brief
    The companion operation to storeParmIndVal.  Takes the array given by that
    function and turns it into the correct return for the current parameter set.
    (ala getVal())
   */
  virtual void getParmDepVal(Wavefunction_data *, Sample_point *,
                             int, //!< electron number
                             Array1 <doublevar> &,
                             //!< values from storeParmIndVal
                             Wf_return &
                             //!< return values, same as from getVal()
                            )=0;

  /*! \brief
    Plots 1d functions from inside the wave function, e.g. constituents
    of the Jastrow factor. Called from PLOT1D method.
  */
  virtual void plot1DInternals(Array1 <doublevar> &,
			       //!< grid of x-values we plot at
			       vector <Array1 <doublevar> > &,
			       //!< one column for each plotted function
			       vector <string> &,
			       //!< description/identification of the column
			       string 
			       //!< prefix for description (BCS, Jastrow, ...)
			       ) {
  }


  /*!
    \brief
    An option for developers to muck around in the internals.


    This isn't really a good function to use, but if you are doing
    development, you can define this for the function you're working
    on to do averages, etc.  This can be removed for release, so
    don't leave it hanging around code.
   */
  virtual void developerAccess(Wavefunction_data *, Sample_point *,
                               Array1 <doublevar> &, Array1 <doublevar> &)
  {
    error("developerAccess isn't defined for this wavefunction."
          "It shouldn't be used unless you're actively working on "
          "something.");
  }

};

int deallocate(Wavefunction * & wfptr);

#include "Sample_point.h"

/*!
\brief
Handles memory management of a sample point and wavefunction
checkpoint.

This helps avoid memory leaks and should simplify
code.
 */
class Storage_container
{
public:
  Storage_container()
  {
    wfStore=NULL;
    sampStore=NULL;
    initialized=0;
  }
  ~Storage_container()
  {
    if(wfStore != NULL)
      delete wfStore;
    if(sampStore != NULL)
      delete sampStore;
  }

  void initialize(Sample_point * sample, Wavefunction * wf)
  {
    assert(wf != NULL);
    if(wfStore != NULL) delete wfStore;
    if(sampStore != NULL ) delete sampStore;
    wfStore=NULL; sampStore=NULL;
    sample->generateStorage(sampStore);
    wf->generateStorage(wfStore);
    initialized=1;
  }

  int isInitialized()
  {
    return initialized;
  }

  void saveUpdate(Sample_point * sample, Wavefunction * wf,
                  int e)
  { 

    wf->saveUpdate(sample, e, wfStore);
    
    sample->saveUpdate(e, sampStore);
  }
  void restoreUpdate(Sample_point * sample, Wavefunction * wf,
                     int e)
  {
    sample->restoreUpdate(e, sampStore);
    wf->restoreUpdate(sample, e, wfStore);
  }
  void restoreUpdate(Sample_point * sample, int e)
  {
    sample->restoreUpdate(e, sampStore);
  }  



private:
  int initialized;

  Wavefunction_storage * wfStore;
  Sample_storage * sampStore;
};



#endif //WAVEFUNCTION_H_INCLUDED

//--------------------------------------------------------------------------

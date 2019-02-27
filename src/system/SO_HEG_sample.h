/*
 
Copyright (C) 2007 Lucas K. Wagner, Jindrich Kolorenc

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


#ifndef SO_HEG_SAMPLE_H_INCLUDED
#define SO_HEG_SAMPLE_H_INCLUDED


#include "Sample_point.h"
#include "SO_HEG_system.h"
class Wavefunction;

class Sample_storage;

/*!
 
*/
class SO_HEG_sample : public Sample_point
{
public:

  SO_HEG_sample() { overall_sign=1.0; overall_phase=0.0; } 
  ~SO_HEG_sample()
  {}

  void init(System * sys);

  void randomGuess();

  int electronSize()
  {
    return nelectrons;
  }
  int ionSize()
  {
    return 0;
  }

  int centerSize()
  {
    return 0;
  }
  doublevar getIonCharge(const int ion)
  {
    return 0;
  }

  virtual int getBounds(Array2 <doublevar> & latvec, Array1 <doublevar> & origin) {
    latvec=parent->latVec;
    origin=parent->origin;
    return 1;
  }
  /*!
    Moves electron e
    position should be an array of length at least 3

  */
  void setElectronPos(const int e,const Array1 <doublevar> & position);
  void setElectronOmega(const int e, const doublevar & omega);
  void translateElectron(const int e, const Array1 <doublevar> & trans);
  
  void getElectronPos(const int e, Array1 <doublevar> & R);
  void getElectronOmega(const int e, doublevar & omega);
  
  void getAllElectronPos(Array2 <doublevar> & pos) {
    assert(pos.GetDim(0) >= nelectrons);
    assert(pos.GetDim(1) >= 3);

    pos=elecpos;
  }

  void moveIon(const int ion, const Array1 <doublevar> & r)
  {
    error("moveIon() not implemented");
  }

  void getIonPos(const int ion, Array1 <doublevar> & r)
  {
    error("getIonPos() not implemented");
  }


  void updateEEDist();
  
  // updateEIDist() is called even if there are no ions defined, therefore
  // has to be empty
  void updateEIDist() { }
  void updateECDist()
  {
    error("updateECDist() not implemented");
  }

  void getECDist(const int e, const int cent,
                 Array1 <doublevar> & distance)
  {
    error("getECDist() not implemented");
  }

  void getEIDist(const int e,const int ion, Array1 <doublevar> & distance)
  {
    error("getEIDist() not implemented");
  }
  /*!
  Returns the vector pointing from e1 to e2.
  */
  void getEEDist(const int e1,const int e2, Array1 <doublevar> & distance)
  {
    assert(distance.GetDim(0) >= 5);
    assert( ! elecDistStale(e1));
    assert( e1 < e2 );
    for(int i=0; i< 5; i++)
    {
      distance(i)=pointdist(e1,e2,i);
    }
  }

  void rawOutput(ostream &);
  void rawInput(istream &);

  void generateStorage(Sample_storage *& );
  void saveUpdate(int, Sample_storage * );
  void restoreUpdate(int, Sample_storage *);

  /*!
  Note that this is not bulletproof; it only works when translateElectron() is used,
  so it doesn't work when you start doing complicated things with setElectronPos(),
  etc.  It's enough for the sampling method, and that's about it.
   */
  doublevar overallSign() { return overall_sign; }
  doublevar overallPhase() { return overall_phase; }

private:

  int nelectrons;
  int dimension;
  Array2 <doublevar> elecpos; //electron positions
  Array1 <doublevar> elecomega;
  Array1 <int> elecDistStale;
  Array3 <doublevar> pointdist;  //this is an upper triangular matrix;
  //the lower is currently wasted
  
  Array2 <doublevar> tmplat; // auxiliary array used in Lucas' updateEEDist()
  Array2 <doublevar> elecpos_lc;   // elecpos in lattice coordinates

  doublevar overall_sign;
  doublevar overall_phase;
  // false for complex-valued wavefunctions, i.e., for non-integer k-points
  bool update_overall_sign;


  SO_HEG_system * parent;     //The System that created this object
};



#endif //SO_HEG_SAMPLE_H_INCLUDED
//-------------------------------------------------------------------------

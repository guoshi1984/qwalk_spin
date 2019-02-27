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

#include "Periodic_sample.h"
#include "ulec.h"
#include "Wavefunction.h"
#include "Periodic_system.h"

//----------------------------------------------------------------------

void Periodic_sample::generateStorage(Sample_storage * & store)
{
  int nions=parent->ions.size();
  store=new Sample_storage;
  store->iondist_temp.Resize(nions,5);
  store->pointdist_temp.Resize(nelectrons,5);
  store->pos_temp.Resize(3);
}

void Periodic_sample::saveUpdate(int e, Sample_storage * store)
{
  getElectronPos(e, store->pos_temp);
  int nions=parent->ions.size();
  for(int d=0; d< 5; d++) {
    for(int i=0; i< nions; i++)
    {
      store->iondist_temp(i,d)=iondist(i,e,d);
    }
    for(int i=0; i<e; i++)
    {
      store->pointdist_temp(i,d)=pointdist(i,e,d);
    }
    for(int j=e; j < nelectrons; j++)
    {
      store->pointdist_temp(j,d)=pointdist(e,j,d);
    }
  }
}

void Periodic_sample::restoreUpdate(int e, Sample_storage * store){
  for(int i=0; i<3; i++)
    elecpos(e,i)=store->pos_temp(i);

  int nions=parent->ions.size();
  for(int d=0; d< 5; d++)  {
    for(int i=0; i< nions; i++)
      iondist(i,e,d)=store->iondist_temp(i,d);

    for(int i=0; i<e; i++)
      pointdist(i,e,d)=store->pointdist_temp(i,d);

    for(int j=e; j < nelectrons; j++)
      pointdist(e,j,d)=store->pointdist_temp(j,d);
  }
  cenDistStale(e)=1; //Let's not save the center distances, but that means
                     //we have to recalculate when we reject.
}


//----------------------------------------------------------------------

void Periodic_sample::updateEIDist() {
  //cout << "Periodic::updateEIDist" << endl;
  //Array1 <doublevar> eivec(3);
  for(int e=0; e< nelectrons; e++) {
    if(ionDistStale(e)) {
      ionDistStale(e)=0;
      int nions=parent->ions.size();
      for(int j=0; j<nions; j++)
        iondist(j,e,1)=1e99;

      //Find the closest image
      Array1 <doublevar> tmpvec(3);
      doublevar tmpdis;
      for(int ion=0; ion < nions; ion++) {
        //for(int d=0; d< 3; d++) eivec(d)=elecpos(e,d)-parent->ions.r(d,ion);
        //cout << "****** ion " << ion << endl;

        for(int aa=-1; aa < 2; aa++) {
          for(int bb=-1; bb < 2; bb++) {
            for(int cc=-1; cc < 2; cc++) {
              tmpdis=0;
              for(int d=0; d< 3; d++) {
                tmpvec(d)=elecpos(e,d)-parent->ions.r(d,ion)
                          +aa*parent->latVec(0,d)
                          +bb*parent->latVec(1,d)
                          +cc*parent->latVec(2,d);
              }
              for(int d=0; d< 3; d++) tmpdis+=tmpvec(d)*tmpvec(d);
              if(tmpdis < iondist(ion,e,1)) {
                //cout << "best one " << tmpdis << endl;
                iondist(ion,e,1)=tmpdis;
                for(int d=0; d< 3; d++) {
                  iondist(ion, e,d+2)=tmpvec(d);
                }
              }
            }
          }
        }
      }
      //cout << "tmpdis " << tmpdis << endl;

      for(int j=0; j<nions; j++)
        iondist(j, e,0)=sqrt(iondist(j, e,1));
    }
  }

  //cout << "done " << endl;

}


//----------------------------------------------------------------------

/*!

 */
void Periodic_sample::updateEEDist() {

  for(int e=0; e< nelectrons; e++) {
    if(elecDistStale(e)==1) {

      elecDistStale(e)=0;
      for(int j=0; j<e; j++) {
        pointdist(j, e,1)=0.0;
      }
      for(int k=e+1; k<nelectrons; k++) {
        pointdist(e,k,1)=0.0;
      }

      //-------------------------
      //Find the closest image
      Array1 <doublevar> tmpvec(3);
      doublevar tmpdis;
      doublevar height2=parent->smallestheight*parent->smallestheight*.25;
      Array2 <doublevar> tmplat(26,3);
      int counter=0;
      for(int aa=-1; aa <= 1; aa++) {
        for(int bb=-1; bb <= 1; bb++) {
          for(int cc=-1; cc <= 1; cc++) {    
            if(aa!=0 || bb !=0 || cc!=0) {
              for(int d=0; d< 3; d++) {
                tmplat(counter,d)=aa*parent->latVec(0,d)
                         +bb*parent->latVec(1,d)
                          +cc*parent->latVec(2,d);
              }
              counter++;
            }
          }
        }
      }

      
      for(int j=0; j< e; j++) {
        //first try within the primitive cell
        for(int d=0; d< 3; d++) 
          pointdist(j,e,d+2)=elecpos(j,d)-elecpos(e,d);
        for(int d=0; d< 3; d++)
          pointdist(j,e,1)+=pointdist(j,e,d+2)*pointdist(j,e,d+2);
        
        //then check outside if we haven't found the minimum image
        if(pointdist(j,e,1) > height2) {
          for(int a=0; a< 26; a++) {
            for(int d=0; d< 3; d++) 
              tmpvec(d)=elecpos(j,d)-elecpos(e,d)+tmplat(a,d);
              
            tmpdis=0;
            for(int d=0; d<3; d++) tmpdis+=tmpvec(d)*tmpvec(d);
            
            if(tmpdis < pointdist(j,e,1)) {
              pointdist(j,e,1)=tmpdis;
              for(int d=0; d< 3; d++) {
                pointdist(j, e, d+2)=tmpvec(d);
              }
            }
            if(tmpdis < height2)
              break;
          }
        }
      }
      
        
       for(int k=e+1; k< nelectrons; k++) {
        //first try within the primitive cell
        for(int d=0; d< 3; d++) 
          pointdist(e,k,d+2)=elecpos(e,d)-elecpos(k,d);
        for(int d=0; d< 3; d++)
          pointdist(e,k,1)+=pointdist(e,k,d+2)*pointdist(e,k,d+2);
        
        //then check outside if we haven't found the minimum image
        if(pointdist(e,k,1) > height2) {
          for(int a=0; a< 26; a++) {
            for(int d=0; d< 3; d++) 
              tmpvec(d)=elecpos(e,d)-elecpos(k,d)+tmplat(a,d);
              
            tmpdis=0;
            for(int d=0; d<3; d++) tmpdis+=tmpvec(d)*tmpvec(d);
            
            if(tmpdis < pointdist(e,k,1)) {
              pointdist(e,k,1)=tmpdis;
              for(int d=0; d< 3; d++) {
                pointdist(e, k, d+2)=tmpvec(d);
              }
            }
            if(tmpdis < height2)
              break;
          }
        }
      }     


      //---------------------------
      for(int j=0; j<e; j++) {
        pointdist(j, e,0)=sqrt(pointdist(j, e,1));
      }
      for(int k=e+1; k<nelectrons; k++) {
        pointdist(e, k,0)=sqrt(pointdist(e, k,1));
        //cout << e << "   " << k << "   " <<  pointdist(0,e,k) << endl;
      }
    }
  }
  //cout << "done" << endl;
}


//----------------------------------------------------------------------
void Periodic_sample::updateECDist() {
  //cout << "periodic_sample::updateECDist() " << endl;
  for(int e=0; e< nelectrons; e++)  {
    if(cenDistStale(e)) {
      cenDistStale(e)=0;
      int ncenters=parent->centerpos.GetDim(0);
      for(int j=0; j<ncenters; j++)
      {
        cendist(e,j,1)=0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<ncenters; j++)
        {
          cendist(e,j, i+2)=elecpos(e,i)-parent->centerpos(j,i);
          cendist(e,j,1)+=cendist(e,j, i+2) * cendist(e,j, i+2);
        }
      }

      for(int j=0; j<ncenters; j++)
        cendist(e,j,0)=sqrt(cendist(e,j,1));
    }
  }
  //cout << "done" << endl;
}
//----------------------------------------------------------------------

void Periodic_sample::init(System * sys) {
  assert(sys != NULL);
  recast(sys, parent); //assign sys to parent


  int nions=parent->ions.size();
  nelectrons=parent->nelectrons(0)+parent->nelectrons(1);

  elecpos.Resize(nelectrons, 3);
  iondist.Resize(nions,nelectrons,5);
  pointdist.Resize(nelectrons, nelectrons,5);

  int ncenters=parent->centerpos.GetDim(0);
  cendist.Resize(nelectrons,ncenters, 5);


  elecDistStale.Resize(nelectrons);
  ionDistStale.Resize(nelectrons);
  cenDistStale.Resize(nelectrons);
  cenDistStale=1;
  elecDistStale=1;
  ionDistStale=1;

  // update_overall_sign is false for complex-valued wavefunctions,
  // i.e., for non-integer k-points
  update_overall_sign=true;
  for (int d=0; d<3; d++) {
    doublevar intpart;
    modf(parent->kpt(d),&intpart);
    if ( parent->kpt(d) != intpart ) update_overall_sign=false;
  }

  /*
  if ( update_overall_sign ) {
    cout << "Sample_point detects real-valued wave function." << endl;
  } else {
    cout << "Sample_point detects complex-valued wave function." << endl;
  }
  */

}


//----------------------------------------------------------------------

void Periodic_sample::randomGuess()
{
  const doublevar range=2.0;  //range of the cube around atoms to fill
  Array1 <doublevar> trialPos(3);

  int e=0;
  //Try to put them around the ions
  //cout << "nions " << ions.size() << endl;
  Array1 <doublevar> ioncenter(3); ioncenter=0;
  for(int ion=0; ion < parent->ions.size(); ion++) {
    for(int d=0;d < 3; d++) ioncenter(d)+=parent->ions.r(d,ion)/parent->ions.size();
  }
  
  for(int spin=0; spin < 2; spin++) {
    for(int ion=0; ion< parent->ions.size(); ion++)
    {
      //cout << "charge " << ions.charge(ion) << endl;
      int laste=e;

      //putting half electrons, rounded up for the first iteration
      //and rounded down for the second one.

      for(;
          e-laste < int(parent->ions.charge(ion)/2)+(1-spin)*int(parent->ions.charge(ion))%2;
          e++)
      {
        if(e< nelectrons) {
          for(int d=0; d< 3; d++)
          {
            trialPos(d)=parent->ions.r(d,ion)+2*range*(rng.ulec()-.5);
          }
          setElectronPos(e,trialPos);
        }
      }
    }
  }

  //if(e > nelectrons)
  //  error("internal:  electrons overflowed in Periodic_sample::randomGuess");

  //Throw the rest randomly around the center of the ions
  for(;e<nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      trialPos(d)=4*range*(rng.ulec()-.5)+ioncenter(d);
    }
    setElectronPos(e, trialPos);
  }

  ionDistStale=1;
  elecDistStale=1;
  cenDistStale=1;

  if(wfObserver)
  {
    wfObserver->notify(all_electrons_move, 0);
  }

}


//----------------------------------------------------------------------

void Periodic_sample::setElectronPos(const int e,
                                   const Array1 <doublevar> & position)
{

  assert( position.GetDim(0) == 3 );
  Array1 <doublevar> temp(position.GetDim(0));
  temp=position;
  Array1 <int> nshift;
  parent->enforcePbc(temp, nshift);

  for(int i=0; i<3; i++)
  {
    //points.r(i,e) = temp(i);
    elecpos(e,i)=temp(i);
  }
  ionDistStale(e)=1;
  elecDistStale(e)=1;
  cenDistStale(e)=1;

  if(wfObserver)
  {
    wfObserver->notify(electron_move, e);
  }

}


void Periodic_sample::translateElectron(const int e, const Array1 <doublevar> & trans) {
  Array1 <doublevar> temp(trans.GetDim(0));
  
  for(int d=0; d< 3; d++) 
    temp(d)=elecpos(e,d)+trans(d);
  Array1<int> nshift;
  parent->enforcePbc(temp, nshift);
  doublevar kdotr=0;
  for(int d=0; d< 3; d++) 
    kdotr+=parent->kpt(d)*nshift(d);
  if ( update_overall_sign ) overall_sign*=cos(pi*kdotr);
  // the sign of the phase here is critical for correct evaluation of
  // pseudopotentials (and density matrices)
  overall_phase-=pi*kdotr;
  // I suppose the value of overall_phase averages around zero, no
  // "modulo" operation should be needed
  //overall_phase=overall_phase - 2*pi * (int) (overall_phase/pi/2);
  
  //if(fabs(kdotr) > 1e-8) 
  //  cout << "shift " << nshift(0) << "   " 
  //      << nshift(1) << "   " << nshift(2) << endl;
  
  for(int i=0; i<3; i++)
    elecpos(e,i)=temp(i);
  
  ionDistStale(e)=1;
  elecDistStale(e)=1;
  cenDistStale(e)=1;

  if(wfObserver)
    wfObserver->notify(electron_move, e);
       
}

//----------------------------------------------------------------------

/*!
\bug
This can possibly screw up if we have more than one sample_point for a system,
and are moving ions for each.  We shouldn't be keeping more than one sample_point
in memory, eventually, anyway, but for now, be careful.  We shouldn't move this
to sys, because we're eventually going to treat the ions as quantum particles..
The solution is to move the ion positions to the sample_point and generalize the
Periodic_system calculation of energies..  Also keep in mind k-points!!
*/
void Periodic_sample::moveIon(const int ion, const Array1 <doublevar> & r) {
  assert(r.GetDim(0) >= 3);


  Array1 <doublevar> temp(r.GetDim(0));
  temp=r;
  Array1<int> nshift;
  parent->enforcePbc(temp, nshift);

  Array1 <doublevar> oldr(3);
  getIonPos(ion, oldr);
  Array1 <doublevar> deltar(3);
  for(int i=0; i< 3; i++) {
    deltar(i)=temp(i)-oldr(i);
  }

  for(int d=0; d< 3; d++) {
    parent->ions.r(d,ion)=temp(d);
    for(int n=0; n< parent->ncenters_atom(ion); n++) {
      parent->centerpos(parent->equiv_centers(ion, n), d)+=deltar(d);
    }
  }

  //int ncenters=parent->centerpos.GetDim(0);
  //for(int c=0; c< ncenters; c++) {
  //  for(int d=0; d< 3; d++)
  //    cout << parent->centerpos(c,d) << "  ";
  //  cout << endl;
  //}

  parent->ion_ewald=parent->ewaldIon();
  //cout << "ewaldIon_done" << endl;
  ionDistStale=1;
  elecDistStale=1;
  cenDistStale=1;

  //cout << "notify " << endl;
  if(wfObserver)
    wfObserver->notify(ion_move, ion);

  //cout << "done " << endl;
}


//----------------------------------------------------------------------

#include <iomanip>

void Periodic_sample::rawOutput(ostream & os)
{
  os << "Periodic_sample\n";
  os << "numIons  " << parent->ions.size() << endl;
  for(int i=0; i< parent->ions.size(); i++)
  {
    for(int d=0; d< 3; d++)
    {
      os << setw(20) << setprecision(12) << parent->ions.r(d,i);
    }
    os << "\n";
  }
  os << "   numElectrons   " << nelectrons << endl;
  for(int e=0; e< nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      os << setw(20) << setprecision(12) << elecpos(e,d);
    }
    os << "\n";
  }
}


void Periodic_sample::rawInput(istream & is)
{
  string text;
  int nions;
  is >> text;
  if(text != "Periodic_sample")
    error("expected Periodic_sample, got ", text);

  is >> text;
  if(text != "numIons")
    error("expected numIons, got ", text);

  is >> nions;
  //cout << "ions" << endl;
  Array2 <doublevar> old_ionpos(nions,3);
  for(int i=0; i< nions; i++)  {
    for(int d=0; d< 3; d++)    {
      is >> old_ionpos(i,d);
    }
    //cout << endl;
  }

  is >> text;
  if(text != "numElectrons")
    error("expected numElectrons, got ", text);
  //cout << "electrons " << endl;
  is >> nelectrons;
  for(int e=0; e< nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      is >> elecpos(e,d);

    }
  }

  //shift the electronic positions if the ions have moved
  updateEIDist();
  for(int e=0; e< nelectrons; e++) {
    Array1 <doublevar> displace(3,0.0);
    doublevar norm=0.0;
    Array1 <doublevar> dist(5);
    for(int at=0; at < nions; at++) {
      getEIDist(e,at, dist);
      doublevar z=getIonCharge(at);
      doublevar func=z/(dist(1)*dist(1));
      norm+=func;      
      //cout << "func " << func << "  dist " << dist(1) << 
      //  "  z " << z << endl;
      for(int d=0; d< 3; d++) 
	displace(d) += func*(parent->ions.r(d,at)-old_ionpos(at,d));
    }
    //cout << "norm = " << norm << endl;
    if(norm>1e-14) {
      for(int d=0; d< 3; d++) 
        elecpos(e,d)+=displace(d)/norm;
      //cout << "displacing electron " << e << " by "
      //     << displace(0)/norm << "   " << displace(1)/norm
      //     << "   " << displace(2)/norm << endl;
    }
  }
  ionDistStale=1;
  elecDistStale=1;
  cenDistStale=1;
  
}

//-------------------------------------------------------------------------

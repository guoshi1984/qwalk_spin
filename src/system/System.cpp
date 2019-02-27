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

#include "System.h"
#include "Molecular_system.h"
#include "Periodic_system.h"
#include "Ring_system.h"
#include "ulec.h"
#include "SHO_system.h"
#include "HEG_system.h"
#include "SO_HEG_system.h"
#include "qmc_io.h"
#include <iostream>

int allocate(vector <string> & syswords,
             System * & sysptr)
{
  assert(sysptr==NULL);
  if(caseless_eq(syswords[0],"MOLECULE")
     || caseless_eq(syswords[0],"MOLECULAR")) 
    sysptr=new Molecular_system;
  
  else if(caseless_eq(syswords[0],"PERIODIC"))
    sysptr=new Periodic_system;
  else if(caseless_eq(syswords[0],"RING"))
    sysptr=new Ring_system;
  else if(caseless_eq(syswords[0],"SHO"))
    sysptr=new SHO_system;
  else if(caseless_eq(syswords[0],"HEG"))
    sysptr=new HEG_system;
  else if(caseless_eq(syswords[0],"SO_HEG"))
  {  sysptr=new SO_HEG_system;
     sysptr->use_spin_sample=1;
  }
  else
    error("Couldn't understand ", syswords[0], " in system section.");
  
  unsigned int pos=0;
  sysptr->read(syswords, pos);
  return 1;
}

//---------------------------------------------------------------------
void System::calcKinetic(Wavefunction_data * wfdata,
                         Sample_point * sample,
                         Wavefunction * wf,
                         Array1 <doublevar> & lap)
{
  assert(lap.GetDim(0)>= wf->nfunc());
  int nelectrons=sample->electronSize();
  int nwf=wf->nfunc();
  lap=0;
  Wf_return temp(nwf,5);
  
  

  Array1 <doublevar> A(3);
  doublevar A2;
  int magnetic_field;
  Array1 <doublevar> modphase(3);
  for(int w=0; w< nwf; w++)
  {
    for(int e=0; e< nelectrons; e++)
    {
     // cout<<"\n calcKinetic ";
      wf->getLap(wfdata, e, temp);
      lap(w)+=temp.amp(w,4);
      //if no magnetic field A=A2=0

    //  cout<<" lap "<<lap(0);

      magnetic_field=getVectorPotential(e, sample, A, A2);

      for(int i=0;i<3;i++)
      {  A(i)=0; modphase(i)=temp.phase(w,i+1)+A(i);}
  //  cout<< " eletron " << e <<" before phase "<< lap(0);  
      if ( temp.is_complex==1 ) {
	lap(w)+= -( modphase(0)*modphase(0)
		   +modphase(1)*modphase(1)
		    +modphase(2)*modphase(2) );//real part should be equal to lap|psi|/|psi|-(grad(phase)+A)^2
      }
      else{
	if(magnetic_field)
	  lap(w)+= -A2;  //real part should be equal to lap|psi|/|psi| -(A)^2	
      }
    //cout<<" electron " << e << " after phase " << lap(0);
   }
      
    lap(w)*=-0.5;
  }
  //cout<<endl;
  //cout << "laplacian " << lap(0) << endl;
  //cout << "Calculating kinetic energy done \n";
}


//----------------------------------------------------------------------
void System::calcKinetic(Wavefunction_data * wfdata,
                         Sample_point * sample,
                         Wavefunction * wf,
                         Array1 <dcomplex> & lap)
{
  assert(lap.GetDim(0)>= wf->nfunc());
  int nelectrons=sample->electronSize();
  int nwf=wf->nfunc();

  lap=dcomplex(0.0,0.0);
  Array1 <doublevar> A(3);
  doublevar A2;
  int magnetic_field;
  Array1 <doublevar> modphase(3);
  
  Wf_return temp(nwf,5);
  for(int w=0; w< nwf; w++)
  {
    for(int e=0; e< nelectrons; e++)
    {

      wf->getLap(wfdata, e, temp);
      lap(w)+=temp.amp(w,4);

      //if no magnetic field A=A2=0
      magnetic_field=getVectorPotential(e, sample, A, A2);

      for(int i=0;i<3;i++)
	modphase(i)=temp.phase(w,i+1)+A(i);
    
      if ( temp.is_complex==1 ) {
	lap(w)+= -( modphase(0)*modphase(0)
		   +modphase(1)*modphase(1)
		    +modphase(2)*modphase(2) ) //real part should be equal to lap|psi|/|psi|-(grad(phase)+A)^2
	  +I*( 2.0*(temp.amp(w,1)*modphase(0)    
		   +temp.amp(w,2)*modphase(1)
		    +temp.amp(w,3)*modphase(2))
	       +temp.phase(w,4) );      //imaginary part  2*(grad(phase)+A)*grad|psi|/|psi| +lap|phase|+....   +grad|A|=0
	       
	
      }
      else{ //temp.is_complex==0
	if(magnetic_field)
	  lap(w)+= -A2 //real part should be equal to lap|psi|/|psi| -(A)^2	
	    +I*(2.0*(temp.amp(w,1)*A(0)    
		     +temp.amp(w,2)*A(1)
		     +temp.amp(w,3)*A(2)));   //imaginary part  2*A*grad|psi|/|psi|
      }
    }
    lap(w)*=-0.5;
  }
  //cout << "laplacian " << lap(0) << endl;
  //cout << "Calculating kinetic energy done \n";
}

//----------------------------------------------------------------------

void System::calcSpinOrbit(Wavefunction_data*, Sample_point*, Wavefunction*,
                           Array1<doublevar>& ){ }

void System::generatePseudo(vector <vector <string> > & words,
                            Pseudopotential * & pseudo) {
  pseudo=new Pseudopotential;
  pseudo->read(words, this);
}



//----------------------------------------------------------------------

int write_xyz(System * sys, ostream & os) {
  
  vector <string> atomlabels;
  sys->getAtomicLabels(atomlabels);
  os<<atomlabels.size()<<endl;
  os << endl; //need a empty line here!!!
  Array1 <doublevar> pos(3);
  for(unsigned int i=0; i<atomlabels.size(); i++) {
    sys->getIonPos(i, pos);
    os<<atomlabels[i] <<" "<< pos(0)
        <<" "<<pos(1)
        <<" "<< pos(2)<<endl;
  }
  return 1;
}

//----------------------------------------------------------------------

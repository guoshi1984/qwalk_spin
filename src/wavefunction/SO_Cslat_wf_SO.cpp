#include <iostream>
#include "Qmc_std.h"
#include "SO_Cslat_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//public
void SO_Cslat_wf::updateSpinOrbit( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);



  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingSO==1) {
      calcSpinOrbit(sample);
      updateEverythingSO=0;
      electronIsStaleSO=0;
    }
   else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {

        if(electronIsStaleSO(e)) {
          assert(!staticSample);
          updateSpinOrbit(sample, e);
          electronIsStaleSO(e)=0;
        }
      }

    }
  }
  //cout<<"updateLap"<<endl; 
}


void SO_Cslat_wf::updateSpinOrbit(Sample_point * sample, int e)
{
  assert(parent!=NULL);
  int s=parent->spin(e);
  sample->updateEIDist();
  dcomplex ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
 
  //check to make sure the determinant isn't zero
  /*for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(detVal(f, det, s)==dcomplex(0.,0.)) {
        cout << "updateLap::WARNING: determinant zero!" << endl;
        Array1 <doublevar> pos(3);
        sample->getElectronPos(e,pos);
        cout << "e= " << e<<" pos " <<  pos(0) << endl;
        calcLap(sample);
        return;
      }
    }
  } */

  parent->cmolecorb->updateSpinOrbit(sample,e ,s,  updatedMoSO);

  for(int d=0 ; d<2 ;d++)
  { for(int i=0; i< updatedMoSO.GetDim(0) ; i++ )
    { MoSOVal(d,e,i)=updatedMoSO(i,d);
    }
  }

}



//private
void SO_Cslat_wf::calcSpinOrbit(Sample_point * sample)
{  
    for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  
    {
      int s=parent->spin(e);
      sample->updateEIDist();
      
      parent->cmolecorb->updateSpinOrbit(sample,e ,s ,  updatedMoSO);

      for(int d=0 ; d<2 ;d++)
      { for(int i=0; i< updatedMoSO.GetDim(0) ; i++ )
        { MoSOVal(d,e,i)=updatedMoSO(i,d);
        }
      }
    } 
     
    int maxmatsize=max(nelectrons(0),nelectrons(1));

    
}


// should call updateLap first to get det value
void SO_Cslat_wf::getSpinOrbit(Wavefunction_data * wfdata,
                          int e, Array1<dcomplex> & value )
{ 
  assert(value.GetDim(0)>=nfunc_);

  int s=parent->spin(e);
  int opp=parent->opspin(e); 
  for(int f=0; f< nfunc_; f++)
  { dcomplex value1(0.0,0.0);
    dcomplex value0(0.0,0.0);  //operator value divided by wavefunction
   
    for(int det=0; det < ndet; det++) 
    { value0+= parent->cdetwt(det)*detVal(f,det,s)*detVal(f,det,opp);
//      cout<<" value0 "<<value0<<endl;
    }
    for(int det=0; det< ndet; det++)
    { dcomplex temp=0;
      for(int j=0; j< nelectrons(s); j++)
      { temp+=(MoSOVal(0,e, parent->occupation(f,det,s)(j))
              +MoSOVal(1,e, parent->occupation(f,det,s)(j)) )
               *inverse(f,det,s)(parent->rede(e),j); 
        
      }  
  //    cout<<" inverse "<<inverse(0,0,0)(0,0) <<endl;  
    //  cout<<" SO "<<MoSOVal(0,0,0)<<" "<<MoSOVal(1,0,0)<<endl;
      value1+=parent->cdetwt(det)*temp*detVal(f,det,s)*detVal(f,det,opp);
    } 
    //cout<<" wavefunction "<<value1;
    value(f)=value1/value0;
//    cout << " hello " << value1 << " "<< value0 << "  "<<  value(f) << endl; 
  }
   
}


void SO_Cslat_wf::updateSpinProj( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

 
  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingSpin==1) {
     calcSpinProj(sample);
      updateEverythingSpin=0;
      electronIsStaleSpin=0;
    }
   else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleSpin(e)) {
          assert(!staticSample);
          updateSpinProj(sample, e);
          electronIsStaleSpin(e)=0;
        }
      }

    }
  }
  //cout<<"updateLap"<<endl; 
}

void SO_Cslat_wf::calcSpinProj(Sample_point * sample)
{
 for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  
 {
   sample->updateEIDist();
   int s=parent->spin(e);

   //for(int f=0; f< nfunc_; f++) {
   //  for(int det=0; det < ndet; det++) {
   //    if(cabs(detVal(f, det, s))==0) {
   //      cout << "updateVal::WARNING: determinant zero!" << endl;
//        calcLap(sample);
   //      return;
   //    }
   //  }
   //} 

  //dcomplex ratio;

   int maxmatsize=max(nelectrons(0),nelectrons(1));
//  Array2 <dcomplex> modet(4,maxmatsize);

  //update all the mo's that we will be using.
  //call the new update orbitals in SO_MO_matrix
   parent->cmolecorb->updateSpinProjVal(sample, e, s,
                              updatedSpinProjMoVal);
   for(int proj=0; proj<4 ;proj++)
   { for (int i=0; i<updatedSpinProjMoVal.GetDim(0); i++)
     { SpinProjMoVal(proj,e,i)=updatedSpinProjMoVal(i,proj) ;}
   }
 }
}


void SO_Cslat_wf::updateSpinProj(Sample_point * sample,int e)
{

  sample->updateEIDist();
  int s=parent->spin(e);

//  for(int f=0; f< nfunc_; f++) {
//    for(int det=0; det < ndet; det++) {
//      if(cabs(detVal(f, det, s))==0) {
//        cout << "updateVal::WARNING: determinant zero!" << endl;
//        calcLap(sample);
//        return;
//      }
//    }
//  }

  //dcomplex ratio;

//  int maxmatsize=max(nelectrons(0),nelectrons(1));
//  Array2 <dcomplex> modet(4,maxmatsize);

  //update all the mo's that we will be using.
  //call the new update orbitals in SO_MO_matrix
  parent->cmolecorb->updateSpinProjVal(sample, e, s,
                              updatedSpinProjMoVal);
  for(int proj=0; proj<4 ;proj++)
  { for (int i=0; i<updatedSpinProjMoVal.GetDim(0); i++)
    { SpinProjMoVal(proj,e,i)=updatedSpinProjMoVal(i,proj) ;}
  }

}

void SO_Cslat_wf::getSpin(Wavefunction_data * wfdata,
                          int e, Array2<dcomplex> & value )
//return sigma_x Psi(R) and sigma_y Psi(R)

{ dcomplex I(0.0,1.0);
  assert(value.GetDim(0)>=nfunc_);
  assert(value.GetDim(1)>=3);
  int s=parent->spin(e);
  int opp=parent->opspin(e); 
  value=0;
  for(int f=0; f< nfunc_; f++)
  { //dcomplex value(0.0,0.0);
    dcomplex value0(0.0,0.0);  //operator value divided by wavefunction
   
    for(int det=0; det < ndet; det++) 
    { value0+= parent->cdetwt(det)*detVal(f,det,s)*detVal(f,det,opp);
     // cout<<" value0 "<<value0<<endl;
    }

   
     for(int det=0; det< ndet; det++)
     { dcomplex temp0(0.0,0.0);
       dcomplex temp1(0.0,0.0);
       for(int j=0; j< nelectrons(s); j++)
       { 
         temp0+=(SpinProjMoVal(2,e, parent->occupation(f,det,s)(j))
              +SpinProjMoVal(3,e, parent->occupation(f,det,s)(j)) )
               *inverse(f,det,s)(parent->rede(e),j); 
         temp1+=(-I*SpinProjMoVal(2,e,parent->occupation(f,det,s)(j))
                +I*SpinProjMoVal(3,e,parent->occupation(f,det,s)(j)) )
               *inverse(f,det,s)(parent->rede(e),j)  ;
        
       }  
  //    cout<<" inverse "<<inverse(0,0,0)(0,0) <<endl;  
    //  cout<<" SO "<<MoSOVal(0,0,0)<<" "<<MoSOVal(1,0,0)<<endl;
  //     cout<<" 0 "<<SpinProjMoVal(0,0,0) << endl;
  //     cout     <<" 1 "<<SpinProjMoVal(1,0,0)<< endl;
  //     cout    <<" i2 "<< I*SpinProjMoVal(2,0,0)<<endl;
  //     cout   <<" i3 "<<I*SpinProjMoVal(3,0,0)<<endl;  
       value(f,0)+=parent->cdetwt(det)*temp0*detVal(f,det,s)*detVal(f,det,opp);
       value(f,1)+=parent->cdetwt(det)*temp1*detVal(f,det,s)*detVal(f,det,opp);
//       cout << " value "<< value(f,1)<<" "<<value0<<endl;
      
       value(f,0)=value(f,0)/value0;
       value(f,1)=value(f,1)/value0;
      } 
   
  }
   
}

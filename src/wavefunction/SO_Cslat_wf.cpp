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
#include <iostream>
#include "Qmc_std.h"
#include "SO_Cslat_wf.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "Slat_wf_data.h"

//----------------------------------------------------------------------


// generateStorage  done
// init        done
// notify      done
// save for static
// saveUpdate
// RestoreUpdate
// updateVal
// updateLap
// storeParmIndVal
// getParmDepVal


// calval
// updateval
// getval
// 
void SO_Cslat_wf::generateStorage(Wavefunction_storage * & wfstore)
{
  //change
  wfstore=new SO_Cslat_wf_storage;
  SO_Cslat_wf_storage * store;
  recast(wfstore, store);
  // 5 component
  //wavefuction storage saves moVal, detVal, inverse
  store->moVal_temp.Resize (5,   nmo);
  for(int i=0; i< 5; i++)
  { for (int mo=0; mo< nmo; mo++)
    { store->moVal_temp(i,mo).Resize(2); 
    }
  } 
  store->detVal_temp.Resize(nfunc_, ndet, 2);
  store->inverse_temp.Resize(nfunc_, ndet, 2);
  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
      for(int s=0; s<2; s++) {
        store->inverse_temp(i,det,s).Resize(nelectrons(s), nelectrons(s));
        store->inverse_temp(i,det,s)=0;
        store->detVal_temp(i,det,s)=1;
      }
      
    }
  }
}


//----------------------------------------------------------------------

void SO_Cslat_wf::init(Wavefunction_data * wfdata)
{
  Slat_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);
  nfunc_=dataptr->nfunc;
  nmo=dataptr->nmo;
  ndet=dataptr->ndet;
//remove
  nelectrons.Resize(2);
//change
  nelectrons=dataptr->nelectrons;
  //cout<<"nelectrons"<<nelectrons(0);  
  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;
 
  spin.Resize(tote);

  spin=dataptr->spin;
  //Properties and intermediate calculation storage.
  moVal.Resize(5,   tote, nmo);
  MoSOVal.Resize(2 , tote, nmo);
  SpinProjMoVal.Resize(4, tote,nmo); 
  for(int i=0; i< 5; i++)
  {  for( int e=0; e< tote; e++)
    { for( int mo=0; mo< nmo; mo++)
      { moVal(i,e,mo).Resize(2); 
      }
    }
  }
  updatedMoVal.Resize(nmo,  5,  2);
  updatedMoSO.Resize(nmo, 2);
  updatedSpinProjMoVal.Resize(nmo,4);
  detVal.Resize (nfunc_, ndet, 2 ) ;

  inverse.Resize(nfunc_, ndet, 2 ) ;

  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
        
     for (int s=0; s<2; s++) {

        inverse(i,det,s).Resize(nelectrons(s), nelectrons(s));
        inverse(i,det,s)=0;
        for(int e=0; e< nelectrons(s); e++) {
          inverse(i,det,s)(e,e)=1;
          inverse(i,det,s)(e,e)=1;
        }

        detVal(i,det,s)=1;
      }

    }
  }


  electronIsStaleVal.Resize(tote);
  electronIsStaleLap.Resize(tote);
  electronIsStaleSO.Resize(tote);
  electronIsStaleSpin.Resize(tote);

  electronIsStaleVal=0;
  electronIsStaleLap=0;
  electronIsStaleSpin=0;
  electronIsStaleSO=0;
  
  updateEverythingVal=1;
  updateEverythingLap=1;
  updateEverythingSO=1;
  updateEverythingSpin=1;
  sampleAttached=0;
  dataAttached=0;
  staticSample=0;
}

//----------------------------------------------------------------------

/*!
Behavior under staticSample:  if sample_static is set, then
we ignore all electron moves in the main algorithm, and the only way
to update based on a move is by the getParmDepVal function.  The
regular update functions will not work in this case.
*/
void SO_Cslat_wf::notify(change_type change, int num)
{
  switch(change)
  {
  case electron_move:
    if(staticSample==0) {
      electronIsStaleVal(num)=1;
      electronIsStaleLap(num)=1;
      electronIsStaleSO(num)=1;
      electronIsStaleSpin(num)=1;
    }
    break;
  case all_electrons_move:
    if(staticSample==0) {
      updateEverythingVal=1;
      updateEverythingLap=1;
      updateEverythingSO=1;
      updateEverythingSpin=1;
    }
    break;
  case wf_parm_change:  
  case all_wf_parms_change:
    if(parent->optimize_mo  ) {
      updateEverythingVal=1;
      updateEverythingLap=1;
      updateEverythingSO=1;
    }
    break;
  case sample_attach:
    sampleAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    updateEverythingSO=1;
    break;
  case data_attach:
    dataAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    updateEverythingSO=1;
    break;
  case sample_static:
    if(!parent->optimize_mo) {
      save_for_static();
      staticSample=1;
    }
    break;
  case sample_dynamic:
    staticSample=0;
    init(parent);
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
    updateEverythingSO=1;
    updateEverythingSpin=1;
  }
}

//----------------------------------------------------------------------

void SO_Cslat_wf::save_for_static() {
  assert(staticSample==0);

  if(!parent->optimize_mo && !parent->optimize_det) {

        
    int totelectrons=nelectrons(0)+nelectrons(1);
    saved_laplacian.Resize(totelectrons, nfunc_, 5);
    Wf_return templap(nfunc_, 5);
    
    for(int e=0; e< totelectrons; e++) {
      getLap(parent, e, templap);
      for(int f=0; f< nfunc_; f++) {
        saved_laplacian(e,f,0)=templap.amp(f,0);
        saved_laplacian(e,f,5)=templap.phase(f,0);
        
      for(int i=1; i< 5; i++) {
        saved_laplacian(e,f,i)=templap.cvals(f,i);
      }
      }
    }
  }

}

//----------------------------------------------------------------------

void SO_Cslat_wf::saveUpdate(Sample_point * sample, int e,
                         Wavefunction_storage * wfstore) {
  if(staticSample==0) {
    SO_Cslat_wf_storage * store;
    recast(wfstore, store);
    int s=spin(e);
    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det<ndet; det++) {
        store->inverse_temp(f,det,s)=inverse(f,det,s);
        store->detVal_temp(f,det,s)=detVal(f,det,s);
      }
    }


    for(int d=0; d< 5; d++) {
      for(int i=0; i< moVal.GetDim(2); i++) {
        store->moVal_temp(d,i)(0)=moVal(d,e,i)(0);
        store->moVal_temp(d,i)(1)=moVal(d,e,i)(1);
      }
    }
  }
}

//----------------------------------------------------------------------

void SO_Cslat_wf::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore)
{

  if(staticSample==0) {
    SO_Cslat_wf_storage * store;
    recast(wfstore, store);
 
   int s=spin(e);

    for(int j=0; j<5; j++) {
      for(int i=0; i<moVal.GetDim(2); i++) {
        moVal(j,e,i)(0)=store->moVal_temp(j,i)(0);
        moVal(j,e,i)(1)=store->moVal_temp(j,i)(1); 
     }
    }
    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det < ndet; det++) {
        inverse(f,det,s)=store->inverse_temp(f,det,s);
        detVal(f,det,s)=store->detVal_temp(f,det,s);
      }
    }

    electronIsStaleVal(e)=0;
    electronIsStaleLap(e)=0;
  }
}

//----------------------------------------------------------------------
//public
void SO_Cslat_wf::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);
  
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingVal==1) {
      calcVal(sample);
      updateEverythingVal=0;
      electronIsStaleVal=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleVal(e)) {
          updateVal(sample, e);
          electronIsStaleVal(e)=0;
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//public
void SO_Cslat_wf::updateLap( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  
  
  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);
  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingLap==1) {
      calcLap(sample);
      updateEverythingVal=0;
      updateEverythingLap=0;
      electronIsStaleLap=0;
      electronIsStaleVal=0;
    }
    else {  
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
                    
        if(electronIsStaleLap(e)) {
          assert(!staticSample);
          updateLap(sample, e);
          electronIsStaleLap(e)=0;
          electronIsStaleVal(e)=0;
        }
      }

    }
  }
  //cout<<"updateLap"<<endl; 
}

//----------------------------------------------------------------------
//Stores wavefunction parameter independent values to the given array. 
void SO_Cslat_wf::storeParmIndVal(Wavefunction_data * wfdata, 
                               Sample_point * sample,
                              int e, Array1 <doublevar> & vals ) {
  assert(vals.GetDim(0) >=2);
  Wf_return newval(nfunc_,1);
  updateVal(wfdata, sample);
  getVal(wfdata, e, newval);
  vals(0)=newval.amp(0,0);
  vals(1)=newval.phase(0,0);

}

//----------------------------------------------------------------------
//The companion operation to storeParmIndVal.
// Takes the array given by that function and 
// turns it into the correct return for the current parameter set. 
void SO_Cslat_wf::getParmDepVal(Wavefunction_data * wfdata,
                            Sample_point * sample,
                            int e,
                            Array1 <doublevar> & oldval,
                            Wf_return & newval) {

  assert(oldval.GetDim(0) >=2);
  assert(newval.amp.GetDim(1) >= 1);
  assert(newval.amp.GetDim(0) >= nfunc_);
  int counter=0;
  newval.amp(0,0)=oldval(counter++);
  newval.phase(0,0)=oldval(counter++);
}


//------------------------------------------------------------------------


void SO_Cslat_wf::calcVal(Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(sample);

}

//------------------------------------------------------------------------
/*!

*/
//private
void SO_Cslat_wf::updateVal(Sample_point * sample,int e)
{

  sample->updateEIDist();
  int s=parent->spin(e);

  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(cabs(detVal(f, det, s))==0) {
        cout << "updateVal::WARNING: determinant zero!" << endl;
        calcLap(sample);
        return;
      }
    }
  }

  dcomplex ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array1 <dcomplex> modet(maxmatsize);

  //update all the mo's that we will be using.
  //call the new update orbitals in SO_MO_matrix
  parent->cmolecorb->updateVal(sample, e, s,
                              updatedMoVal);


  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
   
 //modet is the sum of the two component in the spinor
 for(int i = 0; i < nelectrons(s); i++) {
        modet(i)=updatedMoVal(parent->occupation(f,det,s)(i),0 ,0 )
                +updatedMoVal(parent->occupation(f,det,s)(i),0 ,1 ) ;
      }
//InverseUpdateColumn=A col of A_new * A row of A^-1_old
      ratio=1./InverseUpdateColumn(inverse(f,det,s),//old inverse
                                   modet,//new det
                                   parent->rede(e),// updated col index
                                   nelectrons(s));//matrix dimension;
      detVal(f,det, s)=ratio*detVal(f,det, s);
    }
  }


  for(int i=0; i< updatedMoVal.GetDim(0); i++)
   { for (int spinor=0; spinor<2 ;spinor++) // spinors
     { 
       moVal(0,e,i)(spinor)=updatedMoVal(i,0,spinor);
     }
   }
}

//---------------------------------------------------



/*
void SO_Cslat_wf::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  //Slat_wf_data * slatdata;
  //recast(wfdata, slatdata);

  if(staticSample==0 || parent->optimize_mo ) {
    if(updateEverythingProj==1) {
      calcSpinProjVal(sample);
      updateEverythingProj=0;
      electronIsStaleProj=0;
    }
    else {
      for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
        if(electronIsStaleProj(e)) {
          updateSpinProjVal(sample, e);
          electronIsStaleProj(e)=0;
        }
      }
    }
  }
}

*/

//-----------------------------------

/*void SO_Cslat_wf::calcSpinProjVal(Sample_point * sample)
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


void SO_Cslat_wf::updateSpinProjVal(Sample_point * sample,int e)
{

  sample->updateEIDist();
  int s=parent->spin(e);

  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      if(cabs(detVal(f, det, s))==0) {
        cout << "updateVal::WARNING: determinant zero!" << endl;
//        calcLap(sample);
        return;
      }
    }
  }

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

//  for(int f=0; f< nfunc_; f++)  {
//    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
     //determinant
   
 //modet is the sum of the two component in the spinor
// for(int i = 0; i < nelectrons(s); i++) {
//        for (int j=0; j< 4; j++)
//        { 
//         modet(j,i)=updatedSpinProjMoVal(parent->occupation(f,det,s)(i),j, ) ;
//        }
// }
//InverseUpdateColumn=A col of A_new * A row of A^-1_old
      // for( int j=0; j<4; j++)
      // {
      //  ratio=1./InverseUpdateColumn(inverse(f,det,s),//old inverse
     //                              modet(j),//new det
       //                            parent->rede(e),// updated col index
         //                          nelectrons(s))//matrix dimension;
        //SpinProjDetVal(f,det, s)(j)=ratio*detVal(f,det, s);
      // }
    //}
  //}


  //for(int i=0; i< updatedMoVal.GetDim(0); i++)
  // { for (int spinor=0; spinor<2 ;spinor++) // spinors
  //   { 
  //     moVal(0,e,i)(spinor)=updatedMoVal(i,0,spinor);
  //   }
  // }
}
*/
//------------------------------------------------------------------------


void SO_Cslat_wf::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{

  assert(val.amp.GetDim(0) >= nfunc_);

  Array2 <dcomplex> vals(nfunc_, 1);
  Array1 <doublevar> phase(nfunc_);


  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    for(int f=0; f< nfunc_; f++) {
      vals(f,0)=saved_laplacian(e,f,0);
      phase(f)=saved_laplacian(e,f,5).real();
    }
  }
  else {
    int s=parent->spin(e);
    int opp=parent->opspin(e);
    for(int f=0; f< nfunc_; f++) {

      Array1 <dcomplex> funcval(6);
      funcval=dcomplex(0.0, 0.0);
      for(int det=0; det < ndet; det++) {
        funcval(0) += parent->cdetwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }

      doublevar amp=cabs(funcval(0));


      vals(f,0)=dcomplex(log(amp),0.0);

      //first component of vals: real is the log of amplitude, imag is zero

      if(fabs(funcval(0).imag()) > 1e-8) {
        phase(f)=atan(funcval(0).imag()/funcval(0).real());
	// JK: atan alone has a period of only pi but phase has 2pi period,
	// following should fix this
	if ( funcval(0).real() < 0.0 ) phase(f)+=pi;
      }
      else {
        phase(f)=-.5*pi*(1-sign(funcval(0).real()));
      }
    }
  }

  val.setVals(vals,phase);
}

//don't need-------------------------------------------------------------------
void SO_Cslat_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){

  Array1 <doublevar> phase(1, 0.0);
  Array2 <dcomplex> vals(1,1,dcomplex(0.0,0.0));
  val.setVals(vals, phase);
} 

//-----------------------------------------------------------------------

int SO_Cslat_wf::getParmDeriv(Wavefunction_data *  wfdata, 
			   Sample_point * sample ,
			   Parm_deriv_return & derivatives){

  
/*  int nparms_full=parent->nparms();
  int nparms_start=derivatives.nparms_start;
  int nparms_end=derivatives.nparms_end;
  int nparms=nparms_end-nparms_start;
  
  derivatives.cgradient.Resize(nparms);
  derivatives.chessian.Resize(nparms, nparms);

  derivatives.is_complex=1;
  
  derivatives.gradient.Resize(0);
  derivatives.hessian.Resize(0, 0);
  
  
  if(parent->optimize_mo) {
    error("optimize_mo not implemented for Cslater");
  }
  else if(parent->optimize_det) {
    dcomplex sum=0;
    for(int det=0; det < ndet; det++) {
	sum+=parent->cdetwt(det)*detVal(0,det,0)*detVal(0,det,1);
    }
    
    if(parent->use_csf){
      int counter=0;
      derivatives.cgradient=dcomplex(0.0,0.0);

      for(int csf=0; csf< parent->ncsf; csf++) {
        if(parent->all_weights){
	  if(!parent->complex_detwt){
	    if(csf >= nparms_start && csf <= nparms_end ){
	      dcomplex ctmp_sum=dcomplex(0.0,0.0);
	      for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
		dcomplex ctmp=detVal(0,counter,0)*detVal(0,counter,1);
		ctmp_sum+=parent->CSF(csf)(j)*ctmp;
		counter++;
	      }
	      derivatives.cgradient(csf-nparms_start)=ctmp_sum/sum;
	    }
	    else{
	      counter+=parent->CSF(csf).GetDim(0)-1;
	    }
	  }
	  else{//parent->complex_detwt
	    if(2*csf >= nparms_start && 2*csf+1 <= nparms_end ){
	      dcomplex ctmp_sum=dcomplex(0.0,0.0);
	      for(int j=1;j<parent->CCSF(csf).GetDim(0);j++){
		dcomplex ctmp=detVal(0,counter,0)*detVal(0,counter,1);
		ctmp_sum+=parent->CCSF(csf)(j)*ctmp;
		counter++;
	      }
	      derivatives.cgradient(2*csf-nparms_start)=ctmp_sum/sum;
	      derivatives.cgradient(2*csf+1-nparms_start)=I*(ctmp_sum/sum);//follows from Cauchy-Riemann relations MB
	    }
	    else {
	      counter+=parent->CCSF(csf).GetDim(0)-1;
	    }
	  }
	}//!all_weights
        else{
	  if(!parent->complex_detwt){
	    if(csf > nparms_start && csf <= nparms_end ){
	      dcomplex ctmp_sum=dcomplex(0.0,0.0);
	      for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
		dcomplex ctmp=detVal(0,counter,0)*detVal(0,counter,1);
		ctmp_sum+=parent->CSF(csf)(j)*ctmp;
		counter++;
	      }
	      derivatives.cgradient(csf-1-nparms_start)=ctmp_sum/sum;
	    }
	    else{
	      counter+=parent->CSF(csf).GetDim(0)-1;
	    }
	  }
	  else{//parent->complex_detwt
	    if(2*(csf-1) >= nparms_start && 2*(csf-1)+1 < nparms_end ){
	      dcomplex ctmp_sum=dcomplex(0.0,0.0);
	      for(int j=1;j<parent->CCSF(csf).GetDim(0);j++){
		dcomplex ctmp=detVal(0,counter,0)*detVal(0,counter,1);
		ctmp_sum+=parent->CCSF(csf)(j)*ctmp;
		counter++;
	      }
	      derivatives.cgradient(2*(csf-1)-nparms_start)=ctmp_sum/sum;
	      derivatives.cgradient(2*(csf-1)+1-nparms_start)=I*(ctmp_sum/sum);;//follows from Cauchy-Riemann relations MB
	    }
	    else{
	      counter+=parent->CCSF(csf).GetDim(0)-1;
	    }
	  }
        }
      }//csf
      derivatives.chessian=dcomplex(0.0,0.0);
      return 1;
    }
    else{
      for(int det=0; det < ndet; det++) { 
	dcomplex ctmp=detVal(0,det,0)*detVal(0,det,1);
	
	if(parent->all_weights){
	  if(!parent->complex_detwt){
	    if(det >= nparms_start && det <= nparms_end )
	      derivatives.cgradient(det-nparms_start)=ctmp/sum;
	  }
	  else{
	    if(2*det >= nparms_start && 2*det+1 <= nparms_end ){
	      derivatives.cgradient(2*det-nparms_start)=ctmp/sum; 
	      derivatives.cgradient(2*det+1-nparms_start)=I*(ctmp/sum);//follows from Cauchy-Riemann relations MB
	    }
	  }
	}
	else{
	  if(!parent->complex_detwt){
	    if(det > nparms_start && det <= nparms_end )
	      derivatives.cgradient(det-1-nparms_start)=ctmp/sum;
	  }
	  else{
	    if(2*(det-1) >= nparms_start && 2*(det-1)+1 < nparms_end ){
	      derivatives.cgradient(2*(det-1)-nparms_start)=ctmp/sum;
	      derivatives.cgradient(2*(det-1)+1-nparms_start)=I*(ctmp/sum);//follows from Cauchy-Riemann relations MB
	    }
	  }
	}
      }
      derivatives.chessian=dcomplex(0.0,0.0);
      return 1;
    }
  }
  else { 
    derivatives.cgradient=dcomplex(0.0,0.0);
    derivatives.chessian=dcomplex(0.0,0.0);
    return 1;
  } */
  return 0;
}


//-----------------------------------------------------------------------

void SO_Cslat_wf::getDensity(Wavefunction_data * wfdata, int e,
                         Array2 <doublevar> & dens)
{
/*
  assert(dens.GetDim(0) >= nfunc_);
  //Slat_wf_data * dataptr;
  //recast(wfdata, dataptr);

  int s=parent->spin(e);

  if(ndet > 1)
  {
    error("Haven't done density for several determinants yet.");
  }
  dens=0;
  int det=0;
  for(int f=0; f< nfunc_; f++)
  {
    for(int j=0; j< nelectrons(s); j++)
    {
      doublevar tmp=cabs(moVal(0,e,parent->occupation(f,det,s)(j)));
      dens(f,0)+=tmp*tmp;
    }
  }
*/
}

//----------------------------------------------------------------------------

//private
void SO_Cslat_wf::calcLap(Sample_point * sample)
{

  for(int e=0; e< nelectrons(0)+nelectrons(1); e++)  {
    int s=parent->spin(e);
    sample->updateEIDist();

    //update all the mo's that we will be using, using the lists made in
    //Slat_wf_data(one for each spin).
    parent->cmolecorb->updateLap(sample, e, s,
                                updatedMoVal);
        
    for(int d=0; d< 5; d++)  {
      for(int i=0; i< updatedMoVal.GetDim(0); i++) {
        for(int spinor=0; spinor<2; spinor++)
        {  moVal(d,e,i)(spinor)=updatedMoVal(i,d, spinor);
        }
      }
    }
  }

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  Array2 <dcomplex> modet(maxmatsize, maxmatsize);

  for(int f=0; f< nfunc_; f++)   {
    for(int det=0; det < ndet; det++ ) {
      for(int s=0; s< 2; s++ ) {


        for(int e=0; e< nelectrons(s); e++) {
          int curre=s*nelectrons(0)+e;
          for(int i=0; i< nelectrons(s); i++) {
            modet(e,i)=moVal(0,curre, parent->occupation(f,det,s)(i))(0)
            +moVal(0,curre,parent->occupation(f,det,s)(i))(1);
          }
        }
        
        detVal(f,det,s)=
          TransposeInverseMatrix(modet,inverse(f,det,s), nelectrons(s));

      }
    }
  }
}

//------------------------------------------------------------------------


/*!
\todo
Add a lazy evaluation to this, so it only evaluates the
derivatives if something has really changed.

\bug
When the system is sufficiently large, the wavefunction
is either too small when we begin sampling, or too large
when we actually equilibrate, so we get under/overflows.
Solution:
Put the logarithms further up the chain, since in some
cases, we overflow on the value

*/
//public
void SO_Cslat_wf::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{

  assert(lap.amp.GetDim(1) >= 5);
  assert(lap.amp.GetDim(0) >= nfunc_);

  Array2 <dcomplex> vals(nfunc_, 5,0.0);
  Array1 <doublevar> phase(nfunc_,0.0);
  int zero_func=0;

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
    
   for(int f=0; f< nfunc_; f++) {
      for(int i=0; i< 5; i++)
        vals(f,i)=saved_laplacian(e,f,i);
      phase(f)=saved_laplacian(e,f,5).real();
    }
  }
  else {
    int s=parent->spin(e);
    int opp=parent->opspin(e);

    for(int f=0; f< nfunc_; f++)
    {

      Array1 <dcomplex> funcval(5);
      funcval=dcomplex(0.0, 0.0);

      for(int det=0; det < ndet; det++) {
        funcval(0) += parent->cdetwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      //cout<<" detVal "<<detVal(f,det,s)<<" "<<detVal(f,det,opp);   
      }
      
      for(int i=1; i< 5; i++) {
        for(int det=0; det < ndet; det++) {
          dcomplex temp=0;
          for(int j=0; j<nelectrons(s); j++) 
          {
            temp+=(moVal(i, e, parent->occupation(f,det,s)(j) )(0)
                 + moVal(i, e, parent->occupation(f,det,s)(j) )(1) )
                 *inverse(f,det,s)(parent->rede(e), j);
           //using col of new slater det times the old inverse
          }

          //Prevent catastrophe with a singular matrix.
          //Shouldn't happen much.
          if(detVal(f,det,s)*detVal(f,det,opp)==dcomplex(0.0, 0.0))
          temp=0;

          funcval(i)+=parent->cdetwt(det)*temp
                         *detVal(f,det, s)*detVal(f,det, opp);
        }
      }


      doublevar amp=cabs(funcval(0));
      //cout << "amp " << amp << endl;
      if(fabs(amp) > 1e-16) {
        vals(f,0)=dcomplex(log(amp),0.0);
        if(fabs(funcval(0).imag()) > 1e-8) {
          phase(f)=atan(funcval(0).imag()/funcval(0).real());
	  // JK: atan alone has a period of only pi but phase has 2pi period,
	  // following should fix this
	  if ( funcval(0).real() < 0.0 ) phase(f)+=pi;
        }
        else {
          phase(f)=-.5*pi*(1-sign(funcval(0).real()));
        }
        for(int i=1; i< 5; i++) {
          vals(f,i)=funcval(i)/funcval(0);

        }
      //  cout<<" lap " << funcval(4) <<" " << funcval(0) <<" "<< funcval(4)/funcval(0)<< endl;
       }
      else {
        zero_func=1;
      }
      //cout<< " lap " << vals(0,4);
     
      //cout<<" funcval "<< funcval(4) << " "<<funcval(0);
         
    } //nfunc_
  } //else

  

  if(zero_func) {
        lap.cvals=0;
        lap.amp=0;
        lap.phase=0;
        for(int w=0; w< nfunc_; w++) 
          lap.amp(w,0)=-1e3;
          
  }
  else 
    lap.setVals(vals,phase);
    //cout<<" vals "<< lap.amp(0,4);
    //cout<< " pgs " << lap.phase(0,1)*lap.phase(0,1)+lap.phase(0,2)*lap.phase(0,2)+lap.phase(0,3)*lap.phase(0,3); 

}

//-------------------------------------------------------------------------

/*!
*/

//private , calculate lap of orbitals and the value of determinant 
void SO_Cslat_wf::updateLap(Sample_point * sample,
                        int e ) {
  
  assert(parent != NULL);

  int s=parent->spin(e);
  sample->updateEIDist();
  dcomplex ratio;

  int maxmatsize=max(nelectrons(0),nelectrons(1));
  static Array1 <dcomplex> modet(maxmatsize);

  //check to make sure the determinant isn't zero
  for(int f=0; f< nfunc_; f++) {
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
  }

  
  //update all the mo's that we will be using.
  parent->cmolecorb->updateLap(sample, e, s, updatedMoVal);
  //cout<<"updatedMoVal"<<updatedMoVal(0,0,0); 
  for(int f=0; f<nfunc_; f++)
  {
    for(int det=0; det< ndet; det++)
    {

      //fill the molecular orbitals for this
      //determinant
      for(int i = 0; i < nelectrons(s); i++) {
      
           
       modet(i)=updatedMoVal(parent->occupation(f,det,s)(i), 0 ,0 )
               +updatedMoVal(parent->occupation(f,det,s)(i), 0 , 1 ) ;


//        modet(i)=updatedMoVal(parent->occupation(f,det,s)(i),0,0)+updatedMoVal(parent->occupation(f,det,s)(i),0,1);
      }

      dcomplex tmpratio=InverseUpdateColumn(inverse(f,det,s),
                                             modet,parent->rede(e),
                                             nelectrons(s));
      if(tmpratio==dcomplex(0.,0.))
        ratio=0;
      else ratio=1./tmpratio;

      detVal(f,det, s)=ratio*detVal(f,det, s);
    }
  }


  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal.GetDim(0); i++)
    {
      moVal(d,e,i)(0)=updatedMoVal(i,d,0);
      moVal(d,e,i)(1)=updatedMoVal(i,d,1);
    }


  //cout << "MO determinant " << endl;
  //for(int i=0; i< 2; i++) {
  //  for(int j=0; j< 2; j++) {
  //    cout << moVal(0,i,j) <<"       " ;
  //  }
  //  cout << endl;
  //}

  //cout << "determinant value " << detVal(0,0,0) << endl;

}

/*
void SO_Cslat_wf::getSpinProjVals(Wavefunction_data * wfdata,
                     int e, Wf_return & spv )
{

  assert(spv.spinproj.GetDim(0) >= 4);

  Array2 <dcomplex> vals(nfunc_, 5,0.0);
  Array1 <doublevar> phase(nfunc_,0.0);
  int zero_func=0;

  if(staticSample==1 && parent->optimize_mo==0 && parent->optimize_det==0) {
//    for(int f=0; f< nfunc_; f++) {
//      for(int i=0; i< 5; i++)
//        vals(f,i)=saved_laplacian(e,f,i);
//      phase(f)=saved_laplacian(e,f,5).real();
//    }
  }
  else {
    int s=parent->spin(e);
    int opp=parent->opspin(e);

    for(int f=0; f< nfunc_; f++)
    {

      Array1 <dcomplex> funcval(5);
      funcval=dcomplex(0.0, 0.0);

      for(int det=0; det < ndet; det++) {
        funcval(0) += parent->cdetwt(det)*detVal(f,det,s)*detVal(f,det,opp);
      }

      for(int i=1; i< 4; i++) {
        for(int det=0; det < ndet; det++) {
          dcomplex temp=0;
          for(int j=0; j<nelectrons(s); j++) {
            temp+=SpinProjMoVal(i, e, parent->occupation(f,det,s) ) 
                 *inverse(f,det,s)(parent->rede(e), j);
           //using col of new slater det times the old inverse
          }

          //Prevent catastrophe with a singular matrix.
          //Shouldn't happen much.
          if(detVal(f,det,s)*detVal(f,det,opp)==dcomplex(0.0, 0.0))
          temp=0;

          funcval(i)+=parent->cdetwt(det)*temp
                         *detVal(f,det, s)*detVal(f,det, opp);
        }
      }

    //  doublevar amp=cabs(funcval(0));
      //cout << "amp " << amp << endl;
     // if(fabs(amp) > 1e-16) {
    //    vals(f,0)=dcomplex(log(amp),0.0);
    //    if(fabs(funcval(0).imag()) > 1e-8) {
    //      phase(f)=atan(funcval(0).imag()/funcval(0).real());
	  // JK: atan alone has a period of only pi but phase has 2pi period,
	  // following should fix this
//	  if ( funcval(0).real() < 0.0 ) phase(f)+=pi;
  //      }
    //    else {
    //      phase(f)=-.5*pi*(1-sign(funcval(0).real()));
    //    }
        for(int i=1; i< 4; i++) {
   //       vals(f,i)=funcval(i)/funcval(0);

        }
     // }
      else {
        zero_func=1;
      }
         
    }
  }

 // if(zero_func) {
 //       lap.cvals=0;
 //       lap.amp=0;
 //       lap.phase=0;
 //       for(int w=0; w< nfunc_; w++) 
 //         lap.amp(w,0)=-1e3;
          
//  }
//  else 
//    lap.setVals(vals,phase);
  
 
}
*/
//-------------------------------------------------------------------------

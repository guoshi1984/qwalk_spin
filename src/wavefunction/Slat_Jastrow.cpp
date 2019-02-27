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

#include "Slat_Jastrow.h"
#include "qmc_io.h"
#include "Qmc_std.h"
#include "Sample_point.h"

void Slat_Jastrow::notify(change_type change , int num )
{
  slater_wf->notify(change, num);
  jastrow_wf->notify(change, num);
}



void Slat_Jastrow::init(Wavefunction_data * wfdata)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  dataptr->slater->generateWavefunction(slater_wf);
  dataptr->jastrow->generateWavefunction(jastrow_wf);
          //jastrow points to Jastrow2_wf_data. After the function is called, jastrow_wf points to Jastrow2_wf
  if(jastrow_wf->nfunc() != slater_wf->nfunc())
  {
    error("You must put the same number of functions in the "
          "functions to be multiplied in SLAT-JASTROW.");
  }
  nfunc_=slater_wf->nfunc();
}


void Slat_Jastrow::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new Slat_Jastrow_storage;
  Slat_Jastrow_storage * store;
  recast(wfstore, store);
  slater_wf->generateStorage(store->slat_store);
  jastrow_wf->generateStorage(store->jast_store);
}


void Slat_Jastrow::saveUpdate(Sample_point * sample, int e, Wavefunction_storage * wfstore)
{
  Slat_Jastrow_storage * store;
  recast(wfstore, store);

  slater_wf->saveUpdate(sample, e, store->slat_store);
  jastrow_wf->saveUpdate(sample, e, store->jast_store);

}

void Slat_Jastrow::restoreUpdate(Sample_point * sample, int e, Wavefunction_storage * wfstore)
{
  //cout << "restoreUpdate\n";
  Slat_Jastrow_storage * store;
  recast(wfstore, store);

  slater_wf->restoreUpdate(sample, e, store->slat_store);
  jastrow_wf->restoreUpdate(sample, e, store->jast_store);
}

void Slat_Jastrow::saveUpdate(Sample_point * sample, int e1, int e2, Wavefunction_storage * wfstore)
{
  Slat_Jastrow_storage * store;
  recast(wfstore, store);

  slater_wf->saveUpdate(sample, e1, e2, store->slat_store);
  jastrow_wf->saveUpdate(sample, e1, e2, store->jast_store);

}

void Slat_Jastrow::restoreUpdate(Sample_point * sample, int e1, int e2, Wavefunction_storage * wfstore)
{
  //cout << "restoreUpdate\n";
  Slat_Jastrow_storage * store;
  recast(wfstore, store);

  slater_wf->restoreUpdate(sample, e1, e2, store->slat_store);
  jastrow_wf->restoreUpdate(sample, e1, e2, store->jast_store);
}



void Slat_Jastrow::getVal(Wavefunction_data * wfdata,
                          int e, Wf_return & val)
{
  assert(val.amp.GetDim(0) >= nfunc_);
  Wf_return slat_val(nfunc_, 2);
  Wf_return jast_val(nfunc_, 2);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  
  slater_wf->getVal(dataptr->slater, e, slat_val);
  jastrow_wf->getVal(dataptr->jastrow, e, jast_val);

  if ( slat_val.is_complex==1 || jast_val.is_complex==1 )
    val.is_complex=1;


  for(int i=0; i< nfunc_; i++) {
    val.phase(i,0)=slat_val.phase(i,0)+jast_val.phase(i,0); 
    val.amp(i,0)=slat_val.amp(i,0)+jast_val.amp(i,0);  //add the logarithm
  }
}


void Slat_Jastrow::getSymmetricVal(Wavefunction_data * wfdata,
			     int e, Wf_return & val){

  Wf_return slat_val(nfunc_, 2);
  Wf_return jast_val(nfunc_, 2);
 
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  
  slater_wf->getSymmetricVal(dataptr->slater, e, slat_val);
  jastrow_wf->getSymmetricVal(dataptr->jastrow, e, jast_val);

  //cout <<"slat_val.amp(0,0) "<<slat_val.amp(0,0)<<"  jast_val.amp(0,0) "<<jast_val.amp(0,0)<<endl;

  for(int i=0; i< nfunc_; i++) {
    val.amp(i,0)=slat_val.amp(i,0)+jast_val.amp(i,0);  //add the logarithm
  }
}



void Slat_Jastrow::getDensity(Wavefunction_data * wfdata,
                              int e,
                              Array2 <doublevar> & val)
{
  assert(val.GetDim(0) >= nfunc_);
  Array2 <doublevar> tempval(nfunc_, 2);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  val=0;
  tempval=0;

  slater_wf->getDensity(dataptr->slater, e, tempval);
  jastrow_wf->getDensity(dataptr->jastrow, e, val);

  for(int i=0; i< tempval.GetDim(0); i++)
  {
    val(i,0)*=tempval(i,0);  //multiply..check to see if this is correct
  }
}


void Slat_Jastrow::storeParmIndVal(Wavefunction_data * wfdata,
                                   Sample_point * sample,
                                   int e, Array1 <doublevar> & vals )
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  Array1 <doublevar> slatval(dataptr->slater->valSize());
  Array1 <doublevar> jastval(dataptr->jastrow->valSize());
  slater_wf->storeParmIndVal(dataptr->slater,sample, e, slatval);
  jastrow_wf->storeParmIndVal(dataptr->jastrow,sample, e, jastval);

  int counter=0;
  for(int i=0; i< slatval.GetDim(0); i++)
  {
    vals(counter)=slatval(i);
    counter++;
  }

  for(int i=0; i< jastval.GetDim(0); i++)
  {
    vals(counter)=jastval(i);
    counter++;
  }
}

void Slat_Jastrow::getParmDepVal(Wavefunction_data * wfdata,
                                 Sample_point * sample,
                                 int e,
                                 Array1 <doublevar> & oldval,
                                 Wf_return & newval)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  Array1 <doublevar> slatval(dataptr->slater->valSize());
  Array1 <doublevar> jastval(dataptr->jastrow->valSize());
  int counter=0;
  for(int i=0; i< slatval.GetDim(0); i++)
  {
    slatval(i)=oldval(counter);
    counter++;
  }

  for(int i=0; i< jastval.GetDim(0); i++)
  {
    jastval(i)=oldval(counter);
    counter++;
  }

  Wf_return tempval(nfunc_, 2);

  slater_wf->getParmDepVal(dataptr->slater,sample,  e,slatval, tempval);
  jastrow_wf->getParmDepVal(dataptr->jastrow,sample,  e,jastval, newval);

  for(int i=0; i< tempval.amp.GetDim(0); i++)
  {
    newval.phase(i,0)+=tempval.phase(i,0);  //add phases
    newval.amp(i,0)+=tempval.amp(i,0);  //add the logarithm
  }

  counter=0;
  for(int i=0; i< slatval.GetDim(0); i++)
  {
    oldval(counter)=slatval(i);
    counter++;
  }

  for(int i=0; i< jastval.GetDim(0); i++) {
    oldval(counter)=jastval(i);
    counter++;
  }

}

int Slat_Jastrow::getParmDeriv(Wavefunction_data *  wfdata, 
			       Sample_point * sample ,
			       Parm_deriv_return & derivatives){

  //cout <<"Slat_Jastrow::getParmDeriv"<<endl;
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
   
  Parm_deriv_return slaterval;
  Parm_deriv_return jastval;
  slaterval.need_hessian=jastval.need_hessian=derivatives.need_hessian;
  
  int nslater=dataptr->slater->nparms();
  int njast=dataptr->jastrow->nparms();
  int nparms=derivatives.nparms_end-derivatives.nparms_start;

  if (derivatives.nparms_start<nslater){
    slaterval.nparms_start=derivatives.nparms_start;
    jastval.nparms_start=0;
    
  }
  else{
    jastval.nparms_start=derivatives.nparms_start-nslater;
    slaterval.nparms_start=nslater;
  }
  
  if (derivatives.nparms_end <= nslater){
    slaterval.nparms_end=derivatives.nparms_end;
    jastval.nparms_end=0;
  }
  else{
    jastval.nparms_end=derivatives.nparms_end-nslater;
    slaterval.nparms_end=nslater;
  }
  
  //cout <<slaterval.nparms_start<<" "<<slaterval.nparms_end<<"  "
  //   <<jastval.nparms_start<<" "<<jastval.nparms_end<<endl;
  

  /*
  //new way with usage of extend_parm_deriv; 
  //works as well, you can decrease nparms to be calculated in slater part
  //but not in jastrow part

  Parm_deriv_return retparm;
  retparm.need_hessian=derivatives.need_hessian;
  retparm.nparms_start=derivatives.nparms_start;
  retparm.nparms_end=derivatives.nparms_end;

  if (slaterval.nparms_end-slaterval.nparms_start){
    slater_wf->getParmDeriv(dataptr->slater,sample, slaterval);
    cout <<"adding slater"<<endl;
    extend_parm_deriv(retparm,slaterval);
  }
  
  if (jastval.nparms_end-jastval.nparms_start){
    if(jastval.nparms_end-jastval.nparms_start!=njast)
      error("jastval.nparms_end-jastval.nparms_start!=njast is not supported");
    jastrow_wf->getParmDeriv(dataptr->jastrow,sample, jastval);
    cout<<"adding jastrow"<<endl;
    extend_parm_deriv(retparm,jastval);
  }
  
  derivatives=retparm;
  */

  //good old way, works even for the case jastval.nparms_end-jastval.nparms_start!=njast
  if(nslater)
    slater_wf->getParmDeriv(dataptr->slater,sample, slaterval);
  if(njast){
    jastrow_wf->getParmDeriv(dataptr->jastrow,sample, jastval);
    if(jastval.is_complex)
      error("in  Slat_Jastrow::getParmDeriv we assume Jastrow is always real valued");
  }

  if(!slaterval.is_complex){
    derivatives.gradient.Resize(nparms);
    derivatives.cgradient.Resize(0);
  }
  else {
    derivatives.is_complex=1;
    derivatives.cgradient.Resize(nparms);
    derivatives.gradient.Resize(0);
  }


  for (int i=0;i<nparms;i++){
    if(derivatives.nparms_start+i<nslater)
      if(!slaterval.is_complex)
	derivatives.gradient(i)=slaterval.gradient(i);
      else
	derivatives.cgradient(i)=slaterval.cgradient(i);
    else
      if(!slaterval.is_complex)
	derivatives.gradient(i)=jastval.gradient(i-nslater+derivatives.nparms_start);
      else
	derivatives.cgradient(i)=jastval.gradient(i-nslater+derivatives.nparms_start);
  }



  
  if(derivatives.need_hessian){
    if(!slaterval.is_complex){
      derivatives.hessian.Resize(nparms,nparms);
      derivatives.chessian.Resize(0,0);
    }
    else{
      derivatives.hessian.Resize(0,0);
      derivatives.chessian.Resize(nparms,nparms);
    }
    for(int i=0;i<nparms;i++){
      for(int j=i;j<nparms;j++){
        if(i+derivatives.nparms_start<nslater)
          if (j+derivatives.nparms_start<nslater)
	    if(!slaterval.is_complex)
	      derivatives.hessian(i,j)=slaterval.hessian(i,j);
	    else
	      derivatives.chessian(i,j)=slaterval.chessian(i,j);
          else
	    if(!slaterval.is_complex)
	      derivatives.hessian(i,j)=slaterval.gradient(i)*jastval.gradient(j-nslater+derivatives.nparms_start);
	    else
	     derivatives.chessian(i,j)=slaterval.cgradient(i)*jastval.gradient(j-nslater+derivatives.nparms_start); 
        else
	  if(!slaterval.is_complex)
	    derivatives.hessian(i,j)=jastval.hessian(i-nslater+derivatives.nparms_start,j-nslater+derivatives.nparms_start);
	  else
	    derivatives.chessian(i,j)=jastval.hessian(i-nslater+derivatives.nparms_start,j-nslater+derivatives.nparms_start);
	if(!slaterval.is_complex)
	  derivatives.hessian(j,i)=derivatives.hessian(i,j);
	else
	  derivatives.chessian(j,i)=derivatives.chessian(i,j);
      }
    }
  }
  
  return 1;
}


void Slat_Jastrow::getForceBias(Wavefunction_data * wfdata, int e,
                                Wf_return & bias)
{


  assert(bias.amp.GetDim(0) >= nfunc_);
  assert(bias.amp.GetDim(1) >= 4);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);

  assert(dataptr != NULL);


  Wf_return slat_bias(nfunc_,4);
  Wf_return jast_bias(nfunc_,4);
  slater_wf->getForceBias(dataptr->slater, e, slat_bias);
  jastrow_wf->getForceBias(dataptr->jastrow, e, jast_bias);

  for(int i=0; i< nfunc_; i++)
  {

    //cout << "value slat " << slat_bias(i,0)*exp(slat_bias(i,1) ) << endl;
    //cout << "value jast " << jast_bias(i,0)*exp(jast_bias(i,1) ) << endl;
    bias.phase(i,0)=slat_bias.phase(i,0)+jast_bias.phase(i,0);
    bias.amp(i,0)=slat_bias.amp(i,0)+jast_bias.amp(i,0);
    //cout << "force bias " << bias(0) << endl;
    for(int d=1; d<4; d++)
    {
      bias.amp(i,d)=slat_bias.amp(i,d)+jast_bias.amp(i,d);
      bias.phase(i,d)=slat_bias.phase(i,d)+jast_bias.phase(i,d);
      bias.cvals(i,d)=slat_bias.cvals(i,d)+jast_bias.cvals(i,d);
    }
  }
}


//----------------------------------------------------------------------

/*
void Slat_Jastrow::getLap(Wavefunction_data * wfdata, int e,
                          Wf_return & lap)
{
  //cout << "getLap\n";
  assert(lap.amp.GetDim(0) >= nfunc_);
  assert(lap.amp.GetDim(1) >= 5);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);

  Wf_return slat_lap(nfunc_,5);
  Wf_return jast_lap(nfunc_,5);

  slater_wf->getLap(dataptr->slater, e, slat_lap);
  jastrow_wf->getLap(dataptr->jastrow, e, jast_lap);

  //cout << "slater "; write_array(cout, slat_lap.amp); cout << endl;
  //cout << "jastrow "; write_array(cout, jast_lap.amp); cout << endl;
  //cout << "slater " << slat_lap.amp(0,1) 
  //     << "  " << slat_lap.amp(0,2) << "   " 
  //     << slat_lap.amp(0,3) << "   " << slat_lap.amp(0,4) 
  //     << " jastrow " << jast_lap.amp(0,4) << endl;

  for(int i=0; i< nfunc_; i++)
  {



    //cout << "value slat " << slat_lap(i,0)*exp(slat_lap(i,1) ) << endl;
    lap.phase(i,0)=slat_lap.phase(i,0)+jast_lap.phase(i,0);
    lap.amp(i,0)=slat_lap.amp(i,0)+jast_lap.amp(i,0);
    doublevar dotproduct=0;
    doublevar imagdot=0;
    dcomplex dot(0.0, 0.0);
    for(int d=1; d<4; d++)
    {
      lap.amp(i,d)=slat_lap.amp(i,d)+jast_lap.amp(i,d);
      lap.phase(i,d)=slat_lap.phase(i,d)+jast_lap.phase(i,d);
      dotproduct+=slat_lap.amp(i,d)*jast_lap.amp(i,d)
        +slat_lap.phase(i,d)*jast_lap.phase(i,d);
      imagdot+=slat_lap.amp(i,d)*jast_lap.phase(i,d)
        +slat_lap.phase(i,d)*jast_lap.amp(i,d);


      lap.cvals(i,d)=slat_lap.cvals(i,d)+jast_lap.cvals(i,d);
      dot+=slat_lap.cvals(i,d)*jast_lap.cvals(i,d);

    }
    //cout << "electron " << e << endl;
    //cout << "slat " << slat_lap(i,5)*(-.5) << " jast " << jast_lap(i,5)*(-.5)
    // << endl;
    //cout << "cross product " << dotproduct << endl;
    lap.amp(i,4)=slat_lap.amp(i,4)+jast_lap.amp(i,4)
                  +2*dotproduct;
    lap.phase(i,4)=slat_lap.phase(i,4)+jast_lap.phase(i,4)
      +2.0*imagdot;

    lap.cvals(i,4)=slat_lap.cvals(i,4)+jast_lap.cvals(i,4)
      +2.0*dot;

    //for(int d=0; d< 5; d++) {
    
    //  cout << "i=" <<i << "  d= " << d <<  endl;
    //  cout << "slat_lap " << slat_lap.amp(i,d) << endl;
    //  cout << "jast_lap " << jast_lap.amp(i,d) << endl;
    //  cout << "combined " << lap.amp(i,d) << endl;
      
    //}


  }

}
*/

// JK: complex wavefunctions still behave weird, let's try to rewrite even
// this one, not touching cvals now, though
// The most important is setting is_complex, but I do differ at other places
// too.
void Slat_Jastrow::getLap(Wavefunction_data * wfdata, int e,
                          Wf_return & lap)
{
  //cout << "getLap\n";
  assert(lap.amp.GetDim(0) >= nfunc_);
  assert(lap.amp.GetDim(1) >= 5);
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);

  Wf_return slat_lap(nfunc_,5);
  Wf_return jast_lap(nfunc_,5);

  slater_wf->getLap(dataptr->slater, e, slat_lap);
  jastrow_wf->getLap(dataptr->jastrow, e, jast_lap);

  if ( slat_lap.is_complex==1 || jast_lap.is_complex==1 )
    lap.is_complex=1;

  for(int i=0; i< nfunc_; i++)
  {

    lap.phase(i,0)=slat_lap.phase(i,0)+jast_lap.phase(i,0);
    lap.amp(i,0)=slat_lap.amp(i,0)+jast_lap.amp(i,0);
    doublevar dotproduct=0;
    dcomplex dot(0.0, 0.0);
    for(int d=1; d<4; d++)
    {
      lap.amp(i,d)=slat_lap.amp(i,d)+jast_lap.amp(i,d);
      lap.phase(i,d)=slat_lap.phase(i,d)+jast_lap.phase(i,d);
      dotproduct+=slat_lap.amp(i,d)*jast_lap.amp(i,d);

      lap.cvals(i,d)=slat_lap.cvals(i,d)+jast_lap.cvals(i,d);
      dot+=slat_lap.cvals(i,d)*jast_lap.cvals(i,d);
    }

    lap.amp(i,4)=slat_lap.amp(i,4)+jast_lap.amp(i,4)
      +2*dotproduct;
    //cout<< " " <<e<<" "<< lap.amp(i,4) << " "<<slat_lap.amp(i,4)<<" "<<jast_lap.amp(i,4)<<" "<<dotproduct<<endl;
    lap.phase(i,4)=slat_lap.phase(i,4)+jast_lap.phase(i,4);

    lap.cvals(i,4)=slat_lap.cvals(i,4)+jast_lap.cvals(i,4)
      +2.0*dot;

  }
}


void Slat_Jastrow::getSpinOrbit(Wavefunction_data* wfdata, int e, Array1<dcomplex>& value)
{ 
  dcomplex I(0.0,1.0);
  assert(value.GetDim(0) >= nfunc_);
  
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  Array1<dcomplex> slater_so;
  Array2<dcomplex> slater_spin;
  Wf_return jast_lap(nfunc_,5);
  slater_so.Resize(nfunc_);
  slater_spin.Resize(nfunc_,3); 
  slater_wf->getSpinOrbit(dataptr->slater,e,slater_so);
  slater_wf->getSpin(dataptr->slater,e,slater_spin);
  jastrow_wf->getLap(dataptr->jastrow, e, jast_lap);
  for(int i=0;i<nfunc_;i++ )
  { 
    value(i)=slater_so(i)+slater_spin(i,0)*(-I)*jast_lap.amp(i,2)
            -slater_spin(i,1)*(-I)*jast_lap.amp(i,1);
    //cout<<e <<" x spin "<<  slater_spin(i,0) << " y spin "<<slater_spin(i,1)<<endl;

  }
} 
//----------------------------------------------------------------------

void Slat_Jastrow::updateVal(Wavefunction_data * wfdata, Sample_point * sample)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateVal(dataptr->slater, sample);
  jastrow_wf->updateVal(dataptr->jastrow, sample);
}

void Slat_Jastrow::updateLap(Wavefunction_data * wfdata, Sample_point * sample)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateLap(dataptr->slater, sample);   //slater_wf: pointer to wavefunction class 
                                                  //slater: pointer to Wavefunction_data class
  jastrow_wf->updateLap(dataptr->jastrow, sample);  //jastrow_wf: pointer to wavefunciton class
}

void Slat_Jastrow::updateSpinOrbit(Wavefunction_data* wfdata, Sample_point* sample)
{ Slat_Jastrow_data* dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateSpinOrbit(dataptr->slater, sample);
  slater_wf->updateSpinProj(dataptr->slater,sample);
  jastrow_wf->updateLap(dataptr->jastrow, sample);
}

void Slat_Jastrow::updateForceBias(Wavefunction_data * wfdata, Sample_point * sample)
{
  Slat_Jastrow_data * dataptr;
  recast(wfdata, dataptr);
  assert(dataptr != NULL);
  slater_wf->updateForceBias(dataptr->slater, sample);
  jastrow_wf->updateForceBias(dataptr->jastrow, sample);
}


//----------------------------------------------------------------------


void Slat_Jastrow::plot1DInternals(Array1 <doublevar> & xdata,
				   vector <Array1 <doublevar> > & data,
				   vector <string> & desc,
				   string desc0) {

  slater_wf->plot1DInternals(xdata,data,desc,desc0);
  jastrow_wf->plot1DInternals(xdata,data,desc,desc0);

}


//----------------------------------------------------------------------

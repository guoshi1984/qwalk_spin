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

#include "SO_Split_sampler.h"
#include "twoD_Split_sampler.h"
#include "Split_sample.h"
#include "qmc_io.h"
#include <iostream>
//----------------------------------------------------------------------
/*
void limDrift(Array1 <doublevar> & drift, doublevar tau, drift_type dtype)
{
  doublevar drift2=0;
  doublevar drifta=0;


  if(dtype==drift_cutoff) {
    const doublevar drmax=4;
    
    for(int d=0; d<3; d++)
        drift2+=drift(d)*drift(d);

    drifta=drmax/max( (doublevar) sqrt(drift2), drmax);
    //cout << "drifta " << drifta*tau << endl;
    for(int d=0; d<3; d++) {
        drift(d)*= tau*drifta;
      }
  }
  else if(dtype==drift_cyrus) {
    
    const doublevar acyrus=0.25;
    
    doublevar tau2=tau*2.;
    
    for(int d=0; d<3; d++)
      {
        drift2+=drift(d)*drift(d);
      }
    
    drift2*=acyrus;
    if(drift2 > 1e-16) 
      drifta=(sqrt(1.+tau2*drift2)-1)/(drift2);
    else drifta=0.0;
    
    //cout << "drifta " << drifta << endl;
    for(int d=0; d<3; d++)
      drift(d)*=drifta;
    
  }
  else {
    error("Unknown drift_type in limDrift");
  }
  
}
*/

//-------------------------------------------------------------------

/*!
returns the exponent of the transition probability for
going from point1 to point2
 */
doublevar twoD_Split_sampler::transition_prob(int point1, int point2,
                                         doublevar timestep, 
                                         drift_type dtype) {
  doublevar prob=0;
  Array1 <doublevar> drift(3);
  //cout << "transition probability" << endl;

  drift=trace(point1).drift;

  limDrift(drift, timestep, dtype);

  for(int d=0; d< 3; d++) {
    
    prob-=(trace(point2).pos(d)-trace(point1).pos(d)-drift(d))
      *(trace(point2).pos(d)-trace(point1).pos(d)-drift(d));
  }
  
  //prob-=(trace(point2).omega-trace(point1).omega-drift(4))
  //     *(trace(point2).omega-trace(point1).omega-drift(4));

  prob/=(2.0*timestep);
  
  return prob;
}


//---------------------------------------------------------------------

/*
doublevar runge_kutta_resamp(Point & p1, Point & p2,
                               doublevar timestep, drift_type dtype,
                          int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep, dtype);
  limDrift(dr2, timestep, dtype);
  Array1 <doublevar> avg(3);
  for(int d=0; d< 3; d++) {
    avg(d)=(dr1(d)+dr2(d))/2.0;
  }
  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    green_forward+=(p2.pos(d)-p1.pos(d)-avg(d))*(p2.pos(d)-p1.pos(d)-avg(d));
    //green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))+avg(d)*avg(d);
  }
  return -green_forward/(2*timestep);
}

*/

/*
doublevar runge_kutta_symm(Point & p1, Point & p2,
                               doublevar timestep, drift_type dtype,
                          int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep, dtype);
  limDrift(dr2, timestep, dtype);
  Array1 <doublevar> avg(3);
  for(int d=0; d< 3; d++) {
    avg(d)=(dr1(d)+dr2(d))/2.0;
  }
  //limDrift(avg,timestep, dtype);
  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    //green_forward+=(p2.pos(d)-p1.pos(d)-avg(d))*(p2.pos(d)-p1.pos(d)-avg(d));
    green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))+avg(d)*avg(d);
  }
  return -green_forward/(2*timestep);
}

*/

/*

doublevar linear_symm(Point & p1, Point & p2,
		      doublevar timestep, drift_type dtype,
		      int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep,dtype);
  limDrift(dr2, timestep, dtype);

  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    
    green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))
      +(p2.pos(d)-p1.pos(d))*(dr2(d)-dr1(d))
      +.5*(dr2(d)*dr2(d)+dr1(d)*dr1(d));

    // Runge-kutta move
    //doublevar tmp=p2.pos(d)-p1.pos(d)-.5*(dr2(d)+dr1(d));
    //green_forward+=tmp*tmp;
    //green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))
    //      +.25*(dr2(d)+dr1(d))*(dr2(d)+dr1(d));
  }
  return -green_forward/(2*timestep);
}
*/
//----------------------------------------------------------------------

void twoD_Split_sampler::read(vector <string> & words) {
  unsigned int pos=0;

  readvalue(words, pos=0, divide_, "DIVIDER");
  if(haskeyword(words, pos=0, "RESTRICT_NODES"))
    restrict_nodes=1;
  
  int depth;
  if(readvalue(words, pos=0, depth, "DEPTH")) 
    setRecursionDepth(depth);

  string drifttype;
  if(readvalue(words, pos=0, drifttype, "DRIFT_TYPE")) {
    if(drifttype=="CYRUS")
      setDriftType(drift_cyrus);
    else if(drifttype=="CUTOFF")
      setDriftType(drift_cutoff);
    else error("Didn't understand DRIFT_TYPE ", drifttype);
  }

  
}


/*doublevar Dynamics_generator::greenFunction(Sample_point * sample, Wavefunction * wf,
                                  Wavefunction_data * wfdata, 
                                  Guiding_function * guidingwf,
                                  int e,
                                  Array1 <doublevar> & newpos, 
                                  doublevar timestep,
                                  Dynamics_info & info, Dynamics_info & oldinfo) {
                                    
  //cout << "auxillary " << endl;
  drift_type dtype=drift_cyrus;
  assert(newpos.GetDim(0) >= 3);
  wf->updateLap(wfdata, sample);
  Point p1; p1.lap.Resize(wf->nfunc(), 5);
  p1.sign=sample->overallSign();
  int ndim=sample->ndim();

  sample->getElectronPos(e,p1.pos);
  wf->getLap(wfdata, e, p1.lap);
  guidingwf->getLap(p1.lap, p1.drift);  
  
  
  for(int d=0; d< ndim; d++) 
    p1.translation(d)=newpos(d)-p1.pos(d);
  
  sample->translateElectron(e,p1.translation);
  wf->updateLap(wfdata, sample);
  Point p2; p2.lap.Resize(wf->nfunc(), 5);
  p2.sign=sample->overallSign();
  p2.pos=newpos;
  wf->getLap(wfdata, e, p2.lap);
  guidingwf->getLap(p2.lap, p2.drift);  
  

  info.green_forward=exp(::transition_prob(p1, p2, timestep, dtype));
  info.green_backward=::transition_prob(p2,p1, timestep, dtype);
  info.symm_gf=exp(linear_symm(p1, p2, timestep, dtype));
  
  //cout << "gf:elec " << e << " new wfval " << p2.lap.amp(0,0) << endl;
  info.diffuse_start=p1.drift;
  limDrift(info.diffuse_start, timestep, dtype);
  for(int d=0; d< ndim; d++) {
    //cout << "secondary drift " << info.drift_pos(d) << endl;
    info.diffuse_start(d)+=p1.pos(d);
  }
  info.diffuse_end=newpos;
  info.orig_pos=p1.pos;
  info.new_pos=newpos;
  info.diffusion_rate=0;
  for(int d=0; d< ndim; d++) 
    info.diffusion_rate+=(info.diffuse_end(d)-info.diffuse_start(d))
                        *(info.diffuse_end(d)-info.diffuse_start(d));
  
  info.acceptance=min(1.0, (exp(info.green_backward)/info.green_forward)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap));


  //trying a better gf
  //info.green_forward*=info.acceptance;///oldinfo.acceptance;
  //--------
  

  info.accepted=0;
  doublevar ratio=guidingwf->getTrialRatio(p1.lap,p2.lap)*p1.sign*p2.sign;
  if(ratio < 0) cout << "crossed node " << endl;
  
  return info.acceptance;
}
*/

//----------------------------------------------------------------------
/*!
From x to y in the trace
 */
doublevar twoD_Split_sampler::get_acceptance(Guiding_function * guidingwf
                                        ) {
  //indent = indent + "  ";

  doublevar ratio=guidingwf->getTrialRatio(trace(1).lap, 
                                           trace(0).lap)*trace(0).sign
        *trace(1).sign;

  if(restrict_nodes && ratio < 0) return 0;

  ratio=ratio*ratio;
  //cout << indent << "*********" << endl;
  //cout << indent << "ratio x " << x << " y " << y << " : "  << ratio << endl;

  //int dir=1;
  //if(y-x < 0) {
  //  dir=-1;
  //}


  doublevar prob_transition=0; //exponent of transition probability
  //We use the current move as well(thus y+1 or x-1)
  //for(int i=x+dir; i != y+dir; i+= dir ) {
    //int dist=abs(x-i);
    //cout << indent << "transition from " << y << " to " << y-dir*dist 
    //     << " and " << x << " to " << x+dir*dist << " using " << dist 
    //     << " (timestep " << timesteps(dist)
    //     <<  endl;

    doublevar num=transition_prob(1,0, timesteps(1), dtype);
    doublevar den=transition_prob(0,1, timesteps(1), dtype);
    prob_transition+=num-den;
    //cout << indent <<  "num " << num << " den " << den << " num-den " << num-den << endl;
  //}
  
  //doublevar reject_prob_den=1; //denominator(forward rejections)
  //doublevar reject_prob_num=1; //numerator(backward rejection)
  //for(int i=x+dir; i != y; i+=dir) {
  //  int dist=abs(x-i);
    //cout << indent << "->acceptance from " << y << " to " << y-dir*dist 
    //     << endl;
  //  reject_prob_num*=1-get_acceptance(guidingwf, y, y-dir*dist);
    //cout << indent << "->over " << x << " to " << x+dir*dist << endl;
  //  reject_prob_den*=1-get_acceptance(guidingwf, x,x+dir*dist);
  //}

  doublevar acc=0;

  //const doublevar tiny=1e-10;
  //if(fabs(reject_prob_den) < tiny) {
  //  if(fabs(reject_prob_num) < tiny) acc=0;
  //  else {
      //If we have a zero on the denominator, we know that
      //alpha(y->x)=1/alpha(x->y) must be zero, so therefore the acceptance is 1
    //  acc=1;
      //cout << indent << "warn:denominator zero : " << reject_prob_den << "   num "
      //     << reject_prob_num << endl;
      //error("problem in split_sample; denominator zero");
   // }
  //}
  //else {
    //cout << indent << "prob_transition " << prob_transition
    //    << "   " << reject_prob_num << "   " <<  reject_prob_den << endl;
    acc=ratio*exp(prob_transition);//*reject_prob_num/reject_prob_den;
  //}
  //cout << indent << "acceptance " << acc << endl;
  //indent.erase(indent.end()-1);
  //indent.erase(indent.end()-1);

  return min(1.0,acc);

}

#include "ulec.h"
#include "Wavefunction_data.h"


//----------------------------------------------------------------------
int twoD_Split_sampler::split_driver(int e,
                                Sample_point * sample,
                                Wavefunction * wf, 
                                Wavefunction_data * wfdata,
                                Guiding_function * guidingwf,
                                int depth,
                                Dynamics_info & info,
                                doublevar & efftimestep 
                                ) 
{

  //cout << "primary " << endl;
  
  // in the body of the function, depth =1
  assert(trace.GetDim(0) >= depth);
  assert(recursion_depth_ <= timesteps.GetDim(0));

  if(depth > recursion_depth_) return 0;

  Array1 <doublevar> c_olddrift(3);
  Array1 <doublevar> c_newdrift(3);
  
  c_olddrift=trace(0).drift;  
  limDrift(c_olddrift, timesteps(depth), dtype);

//  int ndim=sample->ndim();

  //cout << "drift " << c_olddrift(0) << "  " 
  //     << c_olddrift(1) << "  " << c_olddrift(2) << endl;
  //info.drift_pos.Resize(3);
  

  // move the sample point using brownian motion with drift
  for(int d=0; d< 2 ; d++) 
  {
    trace(depth).gauss(d)=rng.gasdev();

    trace(depth).translation(d)=trace(depth).gauss(d)*sqrt(timesteps(depth))
        + c_olddrift(d);
//    trace(depth).pos(d)=trace(0).pos(d)
  //      + trace(depth).translation(d);

  }
  
  trace(depth).gauss(2)=0.0;
  trace(depth).translation(2)=0.0;

  //trace(depth).gauss(3)=rng.gasdev();
  //trace(depth).translation(3)=2*trace(depth).gauss(3)*sqrt(timesteps(depth));
  //trace(depth).translation(3)=trace(depth).gauss(3)*sqrt(spintimestep);
  for(int d=0; d<3; d++)
  {  
    trace(depth).pos(d)=trace(0).pos(d)+trace(depth).translation(d);
  
  }
  
 // cout<<"trace before"<<trace(0).pos(2)<<" after "<<trace(1).pos(2)<<endl;  
  /*trace(depth).omega=trace(0).omega+trace(depth).translation(3);
  
  if (trace(depth).omega >= 2*3.1415926 )
  { int t1=trace(depth).omega/(2*3.1415926);
    trace(depth).omega = trace(depth).omega-2*3.1415926*t1; 
  }
  else if ( trace(depth).omega < 0 )
  { int t1=trace(depth).omega/(2*3.1415926);
    trace(depth).omega = trace(depth).omega+2*3.1415926*(t1+1);
  }
  else {}
  */
  doublevar diffusion_rate=0;
  for(int d=0; d< 3; d++) 
  {  diffusion_rate+=trace(depth).gauss(d)*timesteps(depth)*trace(depth).gauss(d);  }
  
  sample->translateElectron(e, trace(depth).translation);
  trace(depth).sign=sample->overallSign();
  
  if(wfdata->supports(laplacian_update) ) {
    wf->updateLap(wfdata, sample);
    wf->getLap(wfdata, e, trace(depth).lap);
  }
  else {
    wf->updateForceBias(wfdata, sample);
    wf->getForceBias(wfdata, e, trace(depth).lap);
  }
  
  guidingwf->getLap(trace(depth).lap, trace(depth).drift);
  //indent="";
  //cout << "#######################acceptance for " << depth << endl;
  doublevar acc=get_acceptance(guidingwf);
  //cout<<"acc "<<acc<<endl; 
  //doublevar acc=get_acceptance(guidingwf);
  //cout << "acceptance for " << depth << " : " <<  acc << endl;    

  info.green_forward=exp(transition_prob(0,depth,timesteps(depth), dtype));
  //cout << "green_forward " << info.green_forward << endl;
  info.green_backward=exp(transition_prob(depth,0,timesteps(depth),dtype));
 
  info.diffusion_rate=diffusion_rate;
  info.acceptance=acc;
  info.orig_pos=trace(0).pos;
  info.diffuse_start.Resize(3);
  for(int d=0; d< 3; d++)
    info.diffuse_start(d)=trace(0).pos(d)+c_olddrift(d);
  //debug_write(cout,"hello\n"); 
  info.diffuse_end=trace(depth).pos;
  info.new_pos=trace(depth).pos;
  info.gauss=trace(depth).gauss;
  
  //info.symm_gf=exp(linear_symm(trace(0), trace(depth), timesteps(depth), dtype));
  //info.resample_gf=info.symm_gf;
  //trying a better gf
  //info.green_forward*=info.acceptance;
  //--------
  
  if (acc+rng.ulec() > 1.0) {
    info.accepted=1;
    return depth;
  }
  else {
    info.accepted=0;
   
    Array1 <doublevar> rev(3,0.0);
    for(int d=0; d< 3; d++) rev(d)=-trace(depth).translation(d);
    sample->translateElectron(e,rev);

    depth++;
    return split_driver(e,sample, wf, wfdata, guidingwf, 
                        depth, info, efftimestep);
  }
   
  //cout<< " position "<<trace(depth).pos(0) <<" "<<trace(depth).pos(1)
    //<<" "<<trace(depth).pos(2) <<endl;
 
}


//----------------------------------------------------------------------

int twoD_Split_sampler::sample(int e,
                          Sample_point * sample, 
                          Wavefunction * wf, 
                          Wavefunction_data * wfdata,
                          Guiding_function * guidingwf,
                          Dynamics_info & info,
                          doublevar & efftimestep
                         ) 
{
  if(! wfStore.isInitialized())
    wfStore.initialize(sample, wf);
  
  wf->updateLap(wfdata, sample);
  wfStore.saveUpdate(sample, wf, e); //wf->saveUpdate , sample->saveUpdate
  trace.Resize(recursion_depth_+1);
  
  for(int i=0; i < recursion_depth_+1; i++) {
    trace(i).lap.Resize(wf->nfunc(), 5);
  }

  timesteps.Resize(recursion_depth_+1);
  timesteps=efftimestep;

    
  for(int i=2; i< recursion_depth_+1; i++) {
    timesteps(i)=efftimestep/pow(divide_,i-1);
  }

  int depth=0;

  sample->getElectronPos(e,trace(depth).pos);
//  sample->getElectronOmega(e, trace(depth).omega);
  wf->getLap(wfdata, e, trace(depth).lap);
  //cout<< " move ";
  trace(depth).sign=sample->overallSign();
  guidingwf->getLap(trace(depth).lap, trace(depth).drift);
  depth++;

  int acc=split_driver(e, sample, wf, wfdata, guidingwf, depth,  
                      info, efftimestep);

  //  cout<< " position "<<trace(depth).pos(0) <<" "<<trace(depth).pos(1)
    //<<" "<<trace(depth).pos(2) <<endl;
 
  if(acc > 0) {
    acceptances(acc-1)++;
    for(int i=0; i< acc; i++) {
      tries(i)++;
    }
  }
  else {
    for(int i=0; i< recursion_depth_; i++) {
      tries(i)++;
    }
  }

  if(!acc) {
    wfStore.restoreUpdate(sample, wf, e);
  }
  //cout << "-----------split done" << endl;
  return acc;
}


//----------------------------------------------------------------------


void twoD_Split_sampler::showStats(ostream & os) {
  doublevar totacc=0;
  for(int i=0; i< recursion_depth_; i++) {
    totacc+=acceptances(i);
  }

  os << "Total acceptance " << totacc/tries(0) << endl;
  for(int i=0; i< recursion_depth_; i++) {
    os << "accept_level" << i << "   " << acceptances(i)/tries(i) 
       << " tries " << tries(i) <<  endl;
  }
}

//----------------------------------------------------------------------


void twoD_Split_sampler::resetStats() {
  acceptances=0;
  tries=0;
}

//----------------------------------------------------------------------
//######################################################################




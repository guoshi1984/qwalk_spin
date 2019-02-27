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
#include "Properties.h"

#include "qmc_io.h"
#include "Wavefunction_data.h"
#include "ulec.h"

//######################################################################

void allocate(vector<string> & words, System * sys, string & runid,
	 Local_density_accumulator *& denspt) { 
  if(words.size() < 1) error("empty density section");
  if(caseless_eq(words[0], "DENSITY"))
    denspt=new One_particle_density;
  else
    error("Didn't understand density keyword",words[0]);
  denspt->init(words, sys, runid);
}

void allocate(vector<string> & words, System * sys, string & runid,
	 Nonlocal_density_accumulator *& nldenspt) { 
  if(words.size() < 1) error("empty non-local density section");
  //  if(caseless_eq(words[0],"OBDM"))
  //    nldenspt=new OBDM;
  //  else if(caseless_eq(words[0],"TBDM"))
  //    nldenspt=new TBDM;
  //  else
  error("Didn't understand non-local density keyword",words[0]);
  nldenspt->init(words, sys, runid);
}

//######################################################################


void One_particle_density::init(vector<string> & words, System * sys,
                                string & runid) {
  int ndim=sys->ndim();
  norm=50;
  outputfile=runid+".cube";
  nup=sys->nelectrons(0);
  unsigned int pos=0;
  
  readvalue(words, pos=0, outputfile,"OUTPUTFILE");
  
  //the domain is chosen as follows:
  //1. If MIN or MAX is input, set the respective quantity to user's values
  //2. If we have a periodic cell, they are set so that the cell is covered.
  //3. Otherwise we set it to include each ion plus  4 bohrs to capture the tails
  min_.Resize(ndim); max.Resize(ndim);

  
  for(int d=0; d< ndim; d++) { 
    min_(d)=0.0; max(d)=0.0;
  }
  
  int nions=sys->nIons();
  atominfo.Resize(nions, 4);
  Array1 <doublevar> ionpos(3);
  for(int i=0; i< nions; i++) {
    sys->getIonPos(i,ionpos);
    for(int d=0; d< 3; d++) {
      atominfo(i,d+1)=ionpos(d);
      if(ionpos(d) < min_(d)) min_(d) = ionpos(d);
      if(ionpos(d) > max(d)) max(d)=ionpos(d);
    }
    atominfo(i,0)=sys->getIonCharge(i);
  }
  
  for(int d=0; d< 3; d++) {
    min_(d)-=4.0;
    max(d)+=4.0;
  }
  
  Array2 <doublevar> latvec;
  if(sys->getBounds(latvec)) { 
    //assume origin is zero for the moment
    min_=0; max=0;
    for(int d=0; d< ndim; d++) {
      for(int i=0; i< ndim; i++) {
        if(latvec(i,d)>0) max(d)+=latvec(i,d);
        if(latvec(i,d)<0) min_(d)+=latvec(i,d);
      }
    }
  }

  
  vector <string> mintxt;
  if(readsection(words, pos=0, mintxt, "MIN")) {
    if(mintxt.size()!=3) error("MIN must have exactly 3 elements");
    for(int d=0;d  < 3; d++) {
      min_(d)=atof(mintxt[d].c_str());
    }
  }

  vector<string> maxtxt;
  if(readsection(words, pos=0, maxtxt, "MAX")) {
    if(maxtxt.size()!=3) error("MAX must have exactly 3 elements");
    for(int d=0;d  < 3; d++) {
      max(d)=atof(maxtxt[d].c_str());
    }
  }

  //------end min/max stuff

  if(!readvalue(words, pos=0, resolution, "RESOLUTION"))
    resolution=.1;
  

  start_electron=0;
  end_electron=sys->nelectrons(0)+sys->nelectrons(1);
  if(haskeyword(words, pos=0, "UP"))
    end_electron=sys->nelectrons(0);
  
  if(haskeyword(words,pos=0, "DOWN"))
    start_electron=sys->nelectrons(0);

  if(end_electron==start_electron)
    error("In one-particle density, UP and DOWN are incompatible");
  
   npoints.Resize(3);
   npoints=1;
  for(int d=0; d < ndim; d++) {
    npoints(d)=int((max(d)-min_(d))/resolution)+1;
  }

  bin.Resize(npoints(0), npoints(1), npoints(2));
  bin=0;
  nsample=0;

  


  //try to recover the previous run's density
  ifstream is(outputfile.c_str());
  if(is) { 
    string dummy;
    //each processor gets an equal share, since we gather them when
    //we write the file
    is >> dummy >> dummy >> nsample;
    nsample /= (mpi_info.nprocs);
    is.ignore(180, '\n'); //finish first line
    is.ignore(180, '\n');
    doublevar dum;
    is >> dum;
    if(dum != nions) error("different number of ions in previous density file");
    for(int d=0; d< 3; d++) {
      is >> dum;
      if(fabs(dum-min_(d)) > 1e-4 )error("different min in density file");
    }
    for(int d=0; d< 3; d++) {
      is >> dum;
      if(fabs(dum-npoints(d))> 1e-4) 
        error("different number of points in density file");
      for(int i=0; i< d; i++) 
        is >> dum;
      is >> dum;
      if(fabs(dum-resolution) > 1e-4) 
        error("different resolution in density file");
      is.ignore(180, '\n');
    }
    for(int at=0; at< nions; at++) 
      is.ignore(180, '\n');
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          is >> bin(x,y,z);

          bin(x,y,z)*=nsample/(norm);
        }
      }
    }   
    
    is.close();
  }
 
}

//----------------------------------------------------------------------

void One_particle_density::accumulate(Sample_point * sample, double weight) { 

  //int nelectrons=sample->electronSize();
  Array1 <int> place(3);
  Array1 <doublevar> epos(3);
  for(int e=start_electron; e< end_electron; e++) {
    sample->getElectronPos(e,epos);
    int use=1;
    for(int d=0; d< 3; d++) {
      place(d)=int( (epos(d)-min_(d))/resolution+0.5);
      if(place(d) <0 || place(d) >= npoints(d))
        use=0;
    }
    if(use) { 
      nsample+=weight;
      bin(place(0), place(1), place(2))+=weight;
    }
  }
}



//----------------------------------------------------------------------

void One_particle_density::write() { 
#ifdef USE_MPI
  Array3 <doublevar> bin_tmp(npoints(0), npoints(1), npoints(2));
  bin_tmp=0;

  //MPI defaults to at max 4MB send queue(at least in the P4 communicator),
  //so we should split the array up into chunks of about 2MB, which is about
  //260,000 doubles, to be safe.  This can be avoided by setting P4_GLOBMEMSIZE
  //to some large value, but it's really annoying to do that, and there's no
  //huge performance benefit.

  
  int n=npoints(0)*npoints(1)*npoints(2);
  double * ptr=bin.v;
  double * ptr2=bin_tmp.v;
  double * last=bin.v+n;
  int interval=260000;
  
  while(ptr<last) { 
    int region=min( int(last-ptr) , interval);
    MPI_Reduce(ptr, ptr2, region, MPI_DOUBLE, MPI_SUM,
	       0, MPI_Comm_grp);
    ptr+=region;
    ptr2+=region;
  }

  doublevar nsample_tmp=parallel_sum(nsample);
#else
  Array3 <doublevar> & bin_tmp(bin);
  doublevar nsample_tmp=nsample;
#endif
  
  if(mpi_info.node==0) {
    ofstream os(outputfile.c_str());
    os << "QWalk: nsamples " << nsample_tmp << "\n";
    os << "Electron density" << endl;
    
    int nions=atominfo.GetDim(0);
    os << "  " << nions << "   " << min_(0) << "   "
       << min_(1) << "   " << min_(2) << endl;
    os << npoints(0) << "   " << resolution << "  0.0000   0.0000" << endl;
    os << npoints(1) << "   0.0000   " << resolution << "  0.0000" << endl;
    os << npoints(2) << "   0.0000    0.0000    " << resolution<< endl;
    
    for(int at=0; at < nions; at++) {
      os << "   " << atominfo(at,0) << "  0.0000  " << atominfo(at,1)
         << "   " << atominfo(at,2) << "   " << atominfo(at,3) << endl;
    }

    int counter=0;
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          os << norm*bin_tmp(x,y,z)/nsample_tmp << "   ";
          if((counter++)%6==5) os << endl;
        }
      }
    }
    os << endl;
    os.close();
    /*  This isn't really used by anyone, is it?  It takes up a lot of space to write the density twice..
    string outputfile2 = outputfile+".dx";
    os.clear();
    os.open(outputfile2.c_str());
    
    /////////////######### DX FILE WRITING ###########
    
    /////##### DX FILE - Probabilities ##########
    //
    //    os << "  " << nions << "   " << min_(0) << "   "
    //       << min_(1) << "   " << min_(2) << endl;
    
    
    os << "object 1 class gridpositions counts " << npoints(0) << " "
      << npoints(1) << " " << npoints(2) << "\n";
    os << "origin " << min_(0) << " "
      << min_(1) << " " << min_(2) << " " << endl;
    
    os << "delta " << resolution << "   0.0   0.0" << endl;
    os << "delta 0.0   " << resolution << "   0.0" << endl;
    os << "delta 0.0   0.0   " << resolution << endl;
    
    os << endl;
    os << "object 2 class gridconnections counts " << npoints(0) << " "
      << npoints(1) << " " << npoints(2) << "\n";
    os << "attribute \"element type\" string \"cubes\" " << endl;
    os << "attribute \"ref\" string \"positions\" " << endl;
    
    
    os << endl;
    os << "object 3 class array type float rank 0 items " << (npoints(0) * npoints(1) * npoints(2)) <<
      " data follows" << endl;
    os << endl;
    
    
    counter=0;
    for(int x=0; x < npoints(0); x++) {
      for(int y=0; y < npoints(1); y++) {
        for(int z=0; z< npoints(2); z++) {
          os << norm*bin_tmp(x,y,z)/nsample_tmp << "   ";
          if((counter++)%6==5) os << endl;
        }
      }
    }
    os << endl;
    
    os << "#attribute \"dep\" string \"positions\" " << endl;
    os << "object \"regular positions regular connections\" class field" << endl;
    os << "component \"positions\" value 1" << endl;
    os << "component \"connections\" value 2" << endl;
    os << "component \"data\" value 3" << endl;
    os << "end" << endl;
    
    
    os.close();
    */
    //////// ######### ion .dx write #############
    /* This needs to utilize the format for .dx code which specifies
      * the attributes of each point based on a non-uniform grid.
      */
    ///////////// ############## end .dx write ########    
    
  }
        
  
}

//######################################################################

//--------------------------------------------------

Properties_manager::Properties_manager() {
  nwf=1;
  maxchildren=3;
  start_avg_step=0;
  current_block=0;
  autocorr_depth=0;
  max_autocorr_depth=10;
  log_file="";
  //maxhist=0;
  is_complex=0;
}

//--------------------------------------------------

Properties_manager::~Properties_manager() {
}

//--------------------------------------------------

void Properties_manager::read(vector <string> & words, 
                              vector <string> & systxt, 
                              vector <string> & wftxt) {

  if(block_avg.GetDim(0) > 0) 
    error("Must call Properties_manager::read before Properties_manager::setSize()");

  unsigned int pos=0;

  readvalue(words, pos=0, max_autocorr_depth, "MAX_AUTOCORR");

}


//--------------------------------------------------

void Properties_manager::setSize(int nwf_, int nblocks, int nsteps, 
                                 int maxwalkers,
                                 System * sys, 
                                 Wavefunction_data * wfdata) {
  assert(nblocks >= 0 );
  assert(nsteps >= 0);
  assert(maxwalkers >= 0);
  assert(nwf_ > 0);
  nwf=nwf_;

  current_block=0;
  
  block_avg.Resize(nblocks);
  for(int i=0; i< nblocks; i++) { 
    if(is_complex){
      block_avg(i).is_complex=is_complex;
      //cout <<" set each block to be complex"<<endl;
    }
    block_avg(i).setSize(nwf, 0,0);
    block_avg(i).aux_size=1;
  }
  
  autocorr_depth=min(nsteps-1, max_autocorr_depth);
  trace.Resize(nsteps, maxwalkers);
}

//--------------------------------------------------




void Properties_manager::insertPoint(int step, 
                                     int walker, 
                                     const Properties_point & pt) {
  assert(walker < trace.GetDim(1));
  assert(step < trace.GetDim(0));
  //cout << "insertpt " << pt.z_pol(0) << endl;
  trace(step, walker)=pt;
 
  if(step >= 1 && pt.parent < 0) {
    error("Problem in insertPoint; parent not set");
  }
}


//--------------------------------------------------


/*!
find the energy autocorrelation.  In general, we have a directed 
graph, so let's start at the bottom and follow parents up.
*/
void Properties_manager::autocorrelation(Array2 <doublevar> & autocorr,
                                         int depth) {

  using namespace Properties_types;
  assert(depth >= 0);
  int nsteps=trace.GetDim(0);
  int nwalkers=trace.GetDim(1);
  int nwf=trace(0,0).kinetic.GetDim(0);
  autocorr.Resize(nwf, depth);
  autocorr=0;
  
  for(int d=1; d< depth+1; d++) {
    int npts=0;
    for(int step=d+start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {

          //Here we trace back up to the dth ancestor of the walker
          int parent=trace(step, walker).parent;
          for(int d1=1; d1 < d; d1++) {
            //if(trace(step-d1, parent).parent <0) {
            //  cout << "step " << step << " d1 " << d1 << "parent " << parent
            //       << endl;
            //}
            assert(trace(step-d1, parent).parent >= 0);
            
            parent=trace(step-d1, parent).parent;
          }
          
          if(trace(step-d, parent).count) {
            npts++;
            
            
            for(int w=0; w< nwf; w++) {
              doublevar en=trace(step, walker).energy(w);
              doublevar enp=trace(step-d, parent).energy(w);

              //doublevar wt_av=block_avg(current_block).weight(w);
              //doublevar wt=trace(step, walker).weight(w);
              //doublevar pwt=trace(step, walker).weight(w);
              doublevar av=block_avg(current_block).avg(total_energy,w);
              
              autocorr(w,d-1)+=(en-av)*(enp-av);
            }


          }
        }
      }
    }
    
    npts=parallel_sum(npts);
    for(int w=0; w< nwf; w++) {
      autocorr(w,d-1)=parallel_sum(autocorr(w,d-1));
      autocorr(w,d-1)/=npts;
    }
  }
  for(int w=0; w< nwf; w++) {
    for(int d=0; d< depth; d++) {
      autocorr(w,d)/= block_avg(current_block).var(total_energy,w);
    }
  }

}


//--------------------------------------------------

/*!
\todo
Use the recurrence relations from 
http://mathworld.wolfram.com/SampleVarianceComputation.html
to remove the extraneous trace in Properties_manager..
 */

void Properties_manager::endBlock() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf, navg_gen);
  assert(trace(0,0).avgrets.GetDim(0)==nwf);
    
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    for(int i=0; i< navg_gen; i++) { 
      block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
      block_avg(current_block).avgrets(w,i).vals=0;
      block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
    }
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          totweight+=wt;
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar wt=trace(step,walker).weight(w);
          avgkin+=trace(step,walker).kinetic(w)*wt;
	  avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w))*wt;
          for(int i=0; i< navg_gen; i++) { 
            for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
              block_avg(current_block).avgrets(w,i).vals(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
            }
          }
        }
      }
    }
    
    weight_sum*=totpts;
    
    block_avg(current_block).totweight=weight_sum;
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;
    block_avg(current_block).avg(ikinetic,w)=0;
    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    
    for(int i=0; i< navg_gen; i++) { 
      for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
        block_avg(current_block).avgrets(w,i).vals(j)=
        parallel_sum(block_avg(current_block).avgrets(w,i).vals(j))/weight_sum;
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;
    doublevar varikin=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
	  doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar en=kin+pot+nonloc;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(ikinetic,w)=0;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    
  }
  

  
  //cout << "autocorrelation " << endl;
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  //cout << "writing block " << endl;
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  final_avg.blockReduce(block_avg, 0, current_block,1);
}



void Properties_manager::endBlockWithSpin() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf, navg_gen);
  assert(trace(0,0).avgrets.GetDim(0)==nwf);
    
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgso=0;
    doublevar avgen=0;
   //average quantities other than energies 
    for(int i=0; i< navg_gen; i++) { 
      block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
      block_avg(current_block).avgrets(w,i).vals=0;
      block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
    }
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          totweight+=wt;
        }
      }
    }
   
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
    
       //cout<<"\n spin_orbit" <<trace(step,walker).spin_orbit(w)<<endl;
        doublevar wt=trace(step,walker).weight(w);
          avgkin+=trace(step,walker).kinetic(w)*wt;
	  avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          //average spin orbit energy
          avgso+=trace(step,walker).spin_orbit(w)*wt;   
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w)
                  +trace(step, walker).spin_orbit(w))*wt;

          for(int i=0; i< navg_gen; i++) { 
            for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
              block_avg(current_block).avgrets(w,i).vals(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
            }
          }
        }
      }
    }
        
    weight_sum*=totpts;
    
    block_avg(current_block).totweight=weight_sum;
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;
    block_avg(current_block).avg(ikinetic,w)=0;
    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(spin_orbit,w)=parallel_sum(avgso)/weight_sum;

    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    
    for(int i=0; i< navg_gen; i++) { 
      for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
        block_avg(current_block).avgrets(w,i).vals(j)=
        parallel_sum(block_avg(current_block).avgrets(w,i).vals(j))/weight_sum;
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varso=0, varen=0;
    doublevar varweight=0;
    doublevar varikin=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
	  doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar so=trace(step,walker).spin_orbit(w);
          doublevar en=kin+pot+nonloc+so;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varso+=(so-block_avg(current_block).avg(spin_orbit,w))
           *(so-block_avg(current_block).avg(spin_orbit,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(ikinetic,w)=0;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    block_avg(current_block).var(spin_orbit,w)=parallel_sum(varso)/weight_sum; 
  }
  

  
  //cout << "autocorrelation " << endl;
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  //cout << "writing block " << endl;
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  final_avg.blockReduce(block_avg, 0, current_block,1);
}





//--------------------------------------------------

/*!
\todo
Use the recurrence relations from 
http://mathworld.wolfram.com/SampleVarianceComputation.html
to remove the extraneous trace in Properties_manager..
 */

void Properties_manager::endBlock_per_step() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf,navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
     

    for(int i=0; i< navg_gen; i++) { 
      block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
      block_avg(current_block).avgrets(w,i).vals=0;
      block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
    }
    
    Array1 <doublevar> weight_per_step(nsteps-start_avg_step,0.0);
    Array1 <int> npts_per_step(nsteps-start_avg_step,0);
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          
          weight_per_step(step-start_avg_step)+=wt;
          npts_per_step(step-start_avg_step)++;
          //cout <<wt<<weight_per_step(step-start_avg_step)<<npts_per_step(step-start_avg_step)<<endl;
          
          totweight+=wt;
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    
    Array1 <int> totnpts_per_step(nsteps-start_avg_step,0);
    Array1 <doublevar> totweight_per_step(nsteps-start_avg_step,0.0);
    for(int step=start_avg_step; step < nsteps; step++) {
      totnpts_per_step(step-start_avg_step)=parallel_sum(npts_per_step(step-start_avg_step));
      totweight_per_step(step-start_avg_step)=parallel_sum(weight_per_step(step-start_avg_step));
      //cout <<"step"<<step<<" totnpts_per_step "<<totnpts_per_step(step-start_avg_step)<<" totweight_per_step " <<totweight_per_step(step-start_avg_step)<<endl;
    }
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    Array1 < Array1 < Array1 <doublevar> > > avgrets_per_step(nsteps-start_avg_step);
    for(int step=start_avg_step; step < nsteps; step++) {
      if(w==0) { 
        avgrets_per_step(step-start_avg_step).Resize(navg_gen);
        for(int i=0; i< navg_gen; i++){
          avgrets_per_step(step-start_avg_step)(i).Resize( block_avg(current_block).avgrets(w,i).vals.GetDim(0));
          avgrets_per_step(step-start_avg_step)(i)=0.0;
        }
      }
      
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar wt=trace(step,walker).weight(w);
          avgkin+=trace(step,walker).kinetic(w)*wt;
          avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w))*wt;
          if(w==0) { 
            for(int i=0; i< navg_gen; i++) { 
              for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
                //block_avg(current_block).avgrets(i).vals(j)+=wt*trace(step,walker).avgrets(i).vals(j);
                avgrets_per_step(step-start_avg_step)(i)(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
              }
            }
          }
        }
      }
    }
    
    weight_sum*=totpts;
    
    //block_avg(current_block).totweight=weight_sum;
    block_avg(current_block).totweight=0.0;
    for(int step=start_avg_step; step < nsteps; step++) {
      if(totnpts_per_step(step-start_avg_step))
	block_avg(current_block).totweight+=totweight_per_step(step-start_avg_step);
    }
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;
    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    //block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    block_avg(current_block).avg(weight,w)=0.0;
    for(int step=start_avg_step; step < nsteps; step++) {
      if(totnpts_per_step(step-start_avg_step))
	block_avg(current_block).avg(weight,w)+=totweight_per_step(step-start_avg_step)/totnpts_per_step(step-start_avg_step);
    }
    block_avg(current_block).avg(weight,w)/=(nsteps-start_avg_step);
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
          //block_avg(current_block).avgrets(i).vals(j)=
          //parallel_sum(block_avg(current_block).avgrets(i).vals(j))/weight_sum;
	  block_avg(current_block).avgrets(w,i).vals(j)=0.0;
	  for(int step=start_avg_step; step < nsteps; step++) {
	    if(totnpts_per_step(step-start_avg_step))
	      block_avg(current_block).avgrets(w,i).vals(j)+=parallel_sum(avgrets_per_step(step-start_avg_step)(i)(j))/totnpts_per_step(step-start_avg_step);
	  }
	  block_avg(current_block).avgrets(w,i).vals(j)/=(nsteps-start_avg_step);
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
          doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar en=kin+pot+nonloc;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
  }
  
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  //cout << "writing block " << endl;
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  
  //cout << "average " << endl;
  final_avg.blockReduce(block_avg, 0, current_block,1);
  //cout << "done" << endl;
}

//--------------------------------------------------


void Properties_manager::endBlockSHDMC() {
  using namespace Properties_types;

  
  int nwalkers=trace.GetDim(1);
  
  int nsteps=trace.GetDim(0);
  int totpts=0;
  //For the generalized averaging, this can be broken if 
  //the calling program isn't consistent with the ordering and 
  //number of Average_returns.  I can't think of any reason why someone
  //would want to do that other than spite, though.
  int navg_gen=trace(0,0).avgrets.GetDim(1);
  block_avg(current_block).avgrets.Resize(nwf,navg_gen);
  
  for(int w=0; w< nwf; w++) {
    doublevar totweight=0;
    doublevar avgkin=0;
    doublevar avgpot=0;
    doublevar avgnonloc=0;
    doublevar avgen=0;
    dcomplex avgckin=dcomplex(0.0,0.0);
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        block_avg(current_block).avgrets(w,i).vals.Resize(trace(0,0).avgrets(w,i).vals.GetDim(0));      
        block_avg(current_block).avgrets(w,i).vals=0;
        block_avg(current_block).avgrets(w,i).type=trace(0,0).avgrets(w,i).type;
      }
    }
    int npts=0;
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          npts++;
          doublevar wt=trace(step,walker).weight(w);
          totweight+=wt;
        }
      }
    }
    //cout << " totpts " << endl;
    totpts=parallel_sum(npts);
    //cout << "donepsum" << endl;
    doublevar weight_sum=parallel_sum(totweight)/totpts;
    
    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar wt=trace(step,walker).weight(w);
	  
	  dcomplex cwt;
	  if(is_complex){
	    cwt=trace(step,walker).cweight(w);
	    avgckin+=trace(step,walker).ckinetic(w)*wt;
	  }          
	  avgkin+=trace(step,walker).kinetic(w)*wt;
          avgpot+=trace(step, walker).potential(w)*wt;
          avgnonloc+=trace(step, walker).nonlocal(w)*wt;
          avgen+=(trace(step, walker).kinetic(w)
                  +trace(step, walker).potential(w)
                  +trace(step, walker).nonlocal(w))*wt;
          if(w==0) { 
	    //using the averaging with w-1 if avgrets(i).type=="linear_delta_der"
            for(int i=0; i< navg_gen; i++) { 
              for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
                if(j<block_avg(current_block).avgrets(w,i).vals.GetDim(0)-1 && block_avg(current_block).avgrets(w,i).type=="linear_delta_der"){
		  if(!trace(step,walker).is_complex)
		    block_avg(current_block).avgrets(w,i).vals(j)+=(wt-1.0)*trace(step,walker).avgrets(w,i).vals(j);
		  else{
		    assert(trace(step,walker).is_complex==1 && is_complex==1);
		    if(j%2==0){
		    //cout <<"j "<<j<<" w "<<cwt.real()<<", "<<cwt.imag()<<" ders "<<trace(step,walker).avgrets(w,i).vals(j)<<" "<<trace(step,walker).avgrets(w,i).vals(j+1)<<endl;
		      block_avg(current_block).avgrets(w,i).vals(j)+=  (cwt.real()-1.0)*trace(step,walker).avgrets(w,i).vals(j)-  cwt.imag()   *trace(step,walker).avgrets(w,i).vals(j+1);
		      block_avg(current_block).avgrets(w,i).vals(j+1)+= cwt.imag()     *trace(step,walker).avgrets(w,i).vals(j)+ (cwt.real()-1)*trace(step,walker).avgrets(w,i).vals(j+1);
		    }
		  }
		}
		else if(j==block_avg(current_block).avgrets(w,i).vals.GetDim(0)-1 && block_avg(current_block).avgrets(w,i).type=="linear_delta_der"){
		  block_avg(current_block).avgrets(w,i).vals(j)+=trace(step,walker).avgrets(w,i).vals(j);
		}
		else 
		  block_avg(current_block).avgrets(w,i).vals(j)+=wt*trace(step,walker).avgrets(w,i).vals(j);
              }
            }
          }
        }
      }
    }
    
    weight_sum*=totpts;
    
    block_avg(current_block).totweight=weight_sum;
    
    block_avg(current_block).avg(kinetic,w)=parallel_sum(avgkin)/weight_sum;

    if(is_complex){
      block_avg(current_block).avg(ikinetic,w)=parallel_sum(avgckin.imag())/weight_sum;
    }
    else
      block_avg(current_block).avg(ikinetic,w)=0.0;

    block_avg(current_block).avg(potential,w)=parallel_sum(avgpot)/weight_sum;
    block_avg(current_block).avg(nonlocal,w)=parallel_sum(avgnonloc)/weight_sum;
    block_avg(current_block).avg(total_energy,w)=parallel_sum(avgen)/weight_sum;
    block_avg(current_block).avg(weight,w)=weight_sum/totpts;
    
    if(w==0) { 
      for(int i=0; i< navg_gen; i++) { 
        for(int j=0; j< block_avg(current_block).avgrets(w,i).vals.GetDim(0); j++) { 
          block_avg(current_block).avgrets(w,i).vals(j)=
          parallel_sum(block_avg(current_block).avgrets(w,i).vals(j))/weight_sum;
        }
      }
    }
    //cout << mpi_info.node << ":npoints" << npts << endl;
    

    //Variance
    doublevar varkin=0, varpot=0, varnonloc=0, varen=0;
    doublevar varweight=0;
    doublevar varikin=0;


    for(int step=start_avg_step; step < nsteps; step++) {
      for(int walker=0; walker < nwalkers; walker++) {
        if(trace(step, walker).count) {
          doublevar kin=trace(step, walker).kinetic(w);
	  doublevar ikin;
	  if(is_complex)
	    ikin=trace(step, walker).ckinetic(w).imag();
          doublevar pot=trace(step, walker).potential(w);
          doublevar nonloc=trace(step, walker).nonlocal(w);
          doublevar en=kin+pot+nonloc;

          doublevar wt=trace(step, walker).weight(w);

          varkin+=(kin-block_avg(current_block).avg(kinetic,w))
            *(kin-block_avg(current_block).avg(kinetic,w))*wt;
	  if(is_complex)
	    varikin+=(ikin-block_avg(current_block).avg(ikinetic,w))
	      *(ikin-block_avg(current_block).avg(ikinetic,w))*wt;
	  
          varpot+=(pot-block_avg(current_block).avg(potential,w))
            *(pot-block_avg(current_block).avg(potential,w))*wt;
          varnonloc+=(nonloc-block_avg(current_block).avg(nonlocal,w))
           *(nonloc-block_avg(current_block).avg(nonlocal,w))*wt;
          varen+=(en-block_avg(current_block).avg(total_energy,w))
            *(en-block_avg(current_block).avg(total_energy,w))*wt;

          varweight+=(wt-block_avg(current_block).avg(weight,w))
            *(wt-block_avg(current_block).avg(weight,w));

        }
      }
    }
    block_avg(current_block).var(kinetic,w)= parallel_sum(varkin)/weight_sum;
    if(is_complex)
      block_avg(current_block).var(ikinetic,w)= parallel_sum(varikin)/weight_sum;
    else
      block_avg(current_block).var(ikinetic,w)=0.0;

    block_avg(current_block).var(potential,w)=parallel_sum(varpot)/weight_sum;
    block_avg(current_block).var(nonlocal, w)=parallel_sum(varnonloc)/weight_sum;
    block_avg(current_block).var(total_energy, w)=parallel_sum(varen)/weight_sum;
    block_avg(current_block).var(weight,w)=parallel_sum(varweight)/totpts;
    
    
  }
  
  autocorrelation(block_avg(current_block).autocorr,
                  autocorr_depth);
  
  if(mpi_info.node==0 && log_file != "") {
    ofstream logout(log_file.c_str(), ios::app);
    logout.precision(16);
    logout << "block { " << endl;
    string indent="   ";
    block_avg(current_block).storeToLog(indent, logout, log_label);
    logout << "} " << endl << endl;
    logout.close();
  }
  current_block++;
  for(int step=0; step < nsteps; step++)
    for(int walker=0; walker < nwalkers; walker++) 
      trace(step, walker).reset();
  
  final_avg.blockReduce(block_avg, 0, current_block,1);
}


//--------------------------------------------------

void Properties_manager::initializeLog(Array2 <Average_generator*> & avg_gen) {
  if(mpi_info.node==0 && log_file != ""){ 
    ofstream os(log_file.c_str(), ios::app);
    os << "init { \n";
    string indent="   ";
    os << indent << "label " << log_label << endl;
    os << indent << "nsys " << avg_gen.GetDim(0) <<  " navg " << avg_gen.GetDim(1) << endl;
    for(int i=0; i< avg_gen.GetDim(0); i++) { 
      for(int j=0; j< avg_gen.GetDim(1); j++) { 
        os << indent << "average_generator { \n";
        avg_gen(i,j)->write_init(indent,os);
        os << indent << "}\n";
      }
    }
    os << "}\n";
    os.close();
  }
}



//--------------------------------------------------

void Properties_manager::initializeLog(Array1 <Average_generator*> & avg_gen) {
  if(mpi_info.node==0 && log_file != ""){ 
    ofstream os(log_file.c_str(), ios::app);
    os << "init { \n";
    string indent="   ";
    os << indent << "label " << log_label << endl;
    for(int i=0; i< avg_gen.GetDim(0); i++) { 
      os << indent << "average_generator { \n";
      avg_gen(i)->write_init(indent,os);
      os << indent << "}\n";
    }
    os << "}\n";
    os.close();
  }
}


//--------------------------------------------------

  
  

void Properties_manager::printBlockSummary(ostream & os) {
  assert(current_block > 0);
  block_avg(current_block-1).printBlockSummary(os);
}

void Properties_manager::printBlockSummaryWithSpin(ostream & os)
{ assert(current_block > 0);
  block_avg(current_block-1).printBlockSummaryWithSpin(os);
}

//--------------------------------------------------


void Properties_manager::printSummary(ostream & os, Array1 <Average_generator*> & avg_gen) {
  if(is_complex)
    final_avg.showComplex(is_complex);
  final_avg.showSummary(os,avg_gen);
}

//--------------------------------------------------

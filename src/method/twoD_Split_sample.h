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


#ifndef twoD_SPLIT_SAMPLE_H_INCLUDED
#define twoD_SPLIT_SAMPLE_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction.h"
#include "System.h"
#include "Sample_point.h"
#include "Guiding_function.h"


/*!
\todo
Make this average the acceptance ratio over all processors

 */
class twoD_Split_sampler:public Dynamics_generator {
 public:

  Split_sampler() {
    divide_=1.0;
    recursion_depth_=1;
    restrict_nodes=0;
    dtype=drift_cyrus;
    //drift_2pt=0;

    acceptances.Resize(recursion_depth_);
    tries.Resize(recursion_depth_);
    acceptances=0;
    tries=0;
  }

  void read(vector <string> & words);

  void setDriftType(drift_type dtype_) {
    dtype=dtype_;
  }
  
  int showinfo(string & indent, ostream & os) {
    os << indent << "recursion depth " << recursion_depth_ << endl;
    os << indent << "timestep divider " << divide_ << endl;
    if(restrict_nodes) os << indent << "restricting node crossings" << endl;
    os << indent << "drift type " << dtype << endl;
    return 1;
  }

  void setDivider(doublevar divide) {
    divide_=divide;
  }

  void setRecursionDepth(int depth) {
    recursion_depth_=depth;
    acceptances.Resize(recursion_depth_);
    tries.Resize(recursion_depth_);
    acceptances=0;
    tries=0;
  }

  int getRecursionDepth() {
    return recursion_depth_;
  }

  /*!
    Returns the step that was accepted, or 0 if 
    there was a rejection.
  */
  int sample(int e,
             Sample_point * sample, 
             Wavefunction * wf, 
             Wavefunction_data * wfdata, 
             Guiding_function * guidewf,
             Dynamics_info & info,
             doublevar & efftimestep
             );
  int split_driver(int e,
                   Sample_point * sample,
                   Wavefunction * wf, 
                   Wavefunction_data * wfdata,
                   Guiding_function * guidewf,
                   int depth,
                   Dynamics_info & info,
                   doublevar & efftimestep);

  //virtual doublevar greenFunction(Sample_point * sample, Wavefunction * wf,
  //                   Wavefunction_data * wfdata, Guiding_function * guidewf,
  //                           int e,
  //                           Array1 <doublevar> & newpos, doublevar timestep,
  //                           Dynamics_info & info);

  doublevar get_acceptance(Guiding_function * guidingwf, int x, int y);
  
  void showStats(ostream & os);
  void resetStats();
 private:
 
  doublevar transition_prob(int point1, int point2,
                            doublevar timestep, 
                            drift_type dtype);
  
  drift_type dtype;
  Array1 <Point> trace;
  Array1 <doublevar> timesteps;
  int recursion_depth_;
  doublevar divide_;
  
  Array1 <doublevar> acceptances;
  Array1 <long int > tries;

  string indent; //for debugging..
  Storage_container wfStore;
};


#endif //2D_SPLIT_SAMPLE_H_INCLUDED

//----------------------------------------------------------------------



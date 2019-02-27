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
#include "Properties_point.h"

void Properties_point::setSize(int nwf) {
  kinetic.Resize(nwf);
  potential.Resize(nwf);
  nonlocal.Resize(nwf);
  weight.Resize(nwf);
  spin_orbit.Resize(nwf);
  wf_val.Resize(nwf, 1);
  weight=1;
  count=0;
  ckinetic.Resize(nwf);
  cweight.Resize(nwf);
  cweight=dcomplex(1.0,0.0);
}


//---------------------------------------------------------------------

void Properties_point::mpiSend(int node) {
#ifdef USE_MPI
  int nwf=kinetic.GetDim(0);

  MPI_Send(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(children.v, children.GetDim(0), MPI_INT,
           node, 0,MPI_Comm_grp);
  MPI_Send(&nchildren, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(&parent, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(&count, 1, MPI_INT, node, 0, MPI_Comm_grp);
  MPI_Send(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp);
  MPI_Send(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_Comm_grp);
  MPI_Send(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp);
  MPI_Send(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_Comm_grp);
  MPI_Send(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_Comm_grp);
  MPI_Send(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE, node, 0, MPI_Comm_grp);
  //added by MB
  MPI_Send(&is_complex, 1, MPI_INT, node, 0, MPI_Comm_grp);
  if(is_complex){
    Array1 <doublevar> realkinetic(nwf);
    Array1 <doublevar> imagkinetic(nwf);
    for(int i=0;i<nwf;i++){
      realkinetic(i)=ckinetic(i).real(); 
      imagkinetic=ckinetic(i).imag();
    }
    MPI_Send(realkinetic.v, realkinetic.GetDim(0), MPI_DOUBLE, 
	     node, 0, MPI_Comm_grp);
    
    MPI_Send(imagkinetic.v, imagkinetic.GetDim(0), MPI_DOUBLE, 
	     node, 0, MPI_Comm_grp); 

     Array1 <doublevar> realweight(nwf);
     Array1 <doublevar> imagweight(nwf);
     for(int i=0;i<nwf;i++){
       realweight(i)=cweight(i).real();
       imagweight(i)=cweight(i).imag();
     }
     MPI_Send(realweight.v, realweight.GetDim(0), MPI_DOUBLE, node, 
	      0, MPI_Comm_grp);
     MPI_Send(imagweight.v, imagweight.GetDim(0), MPI_DOUBLE, node, 
	      0, MPI_Comm_grp);
  }

#else
    error("Properties_point::mpi_send: not using MPI,"
          " this is most likely a bug");
#endif
}

//---------------------------------------------------------------------

void Properties_point::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;

  int nwf;
  MPI_Recv(&nwf, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  //cout << mpi_info.node << "  " << nwf << "  " << naux << endl;

  setSize(nwf);
  MPI_Recv(children.v, children.GetDim(0), MPI_INT,
           node, 0,MPI_Comm_grp, &status);
  MPI_Recv(&nchildren, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(&parent, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);

  MPI_Recv(&count, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  MPI_Recv(kinetic.v, kinetic.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp, & status);

  MPI_Recv(potential.v, potential.GetDim(0), MPI_DOUBLE,
           node, 0, MPI_Comm_grp, & status);

  MPI_Recv(nonlocal.v, nonlocal.GetDim(0), MPI_DOUBLE, 
           node, 0, MPI_Comm_grp, & status);
  MPI_Recv(weight.v, weight.GetDim(0), MPI_DOUBLE, node, 
           0, MPI_Comm_grp, & status);
  MPI_Recv(wf_val.amp.v, wf_val.amp.GetDim(0)*wf_val.amp.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_Comm_grp, & status);
  MPI_Recv(wf_val.phase.v, wf_val.phase.GetDim(0)*wf_val.phase.GetDim(1),
           MPI_DOUBLE,node, 0, MPI_Comm_grp, & status);

  //added by MB
  MPI_Recv(&is_complex, 1, MPI_INT, node, 0, MPI_Comm_grp, &status);
  if(is_complex){
    Array1 <doublevar> realkinetic(nwf);
    Array1 <doublevar> imagkinetic(nwf);

    MPI_Recv(realkinetic.v, realkinetic.GetDim(0), MPI_DOUBLE, 
	     node, 0, MPI_Comm_grp, & status);
    MPI_Recv(imagkinetic.v, imagkinetic.GetDim(0), MPI_DOUBLE, 
	     node, 0, MPI_Comm_grp, & status);
    for(int i=0;i<nwf;i++){
      ckinetic(i).real()=realkinetic(i); 
      ckinetic(i).imag()=imagkinetic(i);
    }
    Array1 <doublevar> realweight(nwf);
    Array1 <doublevar> imagweight(nwf);
  
    MPI_Recv(realweight.v, realweight.GetDim(0), MPI_DOUBLE, node, 
	     0, MPI_Comm_grp, & status);

    MPI_Recv(imagweight.v, imagweight.GetDim(0), MPI_DOUBLE, node, 
	     0, MPI_Comm_grp, & status);
  
    for(int i=0;i<nwf;i++){
      cweight(i).real()=realweight(i);
      cweight(i).imag()=imagweight(i);
    }
  }


#else
  
  error("Properties_point::mpiRecieve: not using MPI,"
        " this is most likely a bug");
#endif
}

//----------------------------------------------------------------------

#include "qmc_io.h"

void Properties_point::read(istream & is) { 
  int nwf;
  string dummy;
  const string errmsg="Properties error in checkpoint read";
  is >> dummy >> nwf;
  if(dummy != "nwf") error("expected nwf, got ", dummy);
  is >> dummy;
  if(dummy=="naux") {
    debug_write(cout,"Trying to read from old Properties_point..");
    int naux;
    is >> naux;
    if(naux > 0) error("Don't support auxiliary wave functions any more!");
    is >> dummy >> dummy;
  }
  setSize(nwf);
  is >> dummy;
  read_array(is, nwf, kinetic);
  is >> dummy;
  read_array(is,nwf, potential);
  is >> dummy;
  read_array(is, nwf, nonlocal);
  is >> dummy;
  read_array(is, nwf, weight);
  is >> dummy >> dummy;  //Wf_val { 
  wf_val.read(is);
  is >> dummy; //}
}

//----------------------------------------------------------------------


void Properties_point::write(string & indent, ostream & os) { 
  int nwf=kinetic.GetDim(0);
  os << indent << "nwf " << nwf << endl;
  os << indent << "kinetic ";
  write_array(os, kinetic);
  os << endl << indent << "potential ";
  write_array(os, potential);
  os << endl << indent << "nonlocal ";
  write_array(os, nonlocal);
  os << endl << indent << "weight ";
  write_array(os, weight);
  os << endl;

  os << indent << "Wf_val { \n";
  string indent2=indent+"  ";
  wf_val.write(indent2, os);
  os << indent << "}\n";
}

//----------------------------------------------------------------------

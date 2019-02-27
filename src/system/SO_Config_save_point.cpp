#include "SO_Config_save_point.h"

void SO_Config_save_point::write(ostream & os)
{ os<<"nElec "<<electronpos.GetDim(0)<<" ndim "<<4<<endl;
  for(int e=0; e<electronpos.GetDim(0); e++)
  { for(int d=0; d<3; d++)
    {  os<<electronpos(e)(d)<<" ";
    }
    os<<electronomega(e)<< " ";
    os<<endl;
  }
  os << endl;
}

void SO_Config_save_point::read(istream & is)
{ 
  string dummy;
  is >> dummy;
  assert(dummy=="nElec");
  int nelec;
  is >> nelec;
  electronpos.Resize(nelec);
  electronomega.Resize(nelec);
  is >> dummy;
  assert(dummy=="ndim");
  int ndim;
  is >> ndim;
  for(int e=0; e< nelec; e++) 
  {
    electronpos(e).Resize(3);
    
    for(int d=0; d< 3; d++) 
    {
      is >> electronpos(e)(d);
    }
    is >> electronomega(e);
      
  }

}
 
void SO_Config_save_point::savePos(Sample_point * sample) {

  int nelectrons=sample->electronSize();

  if(electronpos.GetDim(0)!=nelectrons) {
    electronpos.Resize(nelectrons);
    for(int i=0; i< nelectrons; i++) {
      electronpos(i).Resize(3);

    }
  }


  for(int i=0; i< nelectrons; i++) {
    sample->getElectronPos(i,electronpos(i));
  }


}

void SO_Config_save_point::restorePos(Sample_point * sample) 
{ 
  int nelectrons=sample->electronSize();
  assert(nelectrons==electronpos.GetDim(0));
  for(int i=0; i< nelectrons; i++)
    sample->setElectronPos(i,electronpos(i));
}

void SO_Config_save_point::getPos(int e, Array1<doublevar> & r)
{  r=electronpos(e);
}

void SO_Config_save_point::saveOmega(Sample_point * sample)
{ int nelectrons=sample->electronSize();
  if(electronomega.GetDim(0)!=nelectrons )
  { electronomega.Resize(nelectrons);
  }
  
  for(int i=0; i< nelectrons; i++)
  { sample->getElectronOmega(i,electronomega(i)); }
  
}


void SO_Config_save_point::restoreOmega(Sample_point * sample)
{ int nelectrons=sample->electronSize();
  assert(nelectrons==electronomega.GetDim(0));
  for (int i=0; i<nelectrons; i++)
   sample-> setElectronOmega(i,electronomega(i));
}


void SO_Config_save_point::getOmega(int e, doublevar & omega)
{ omega=electronomega(e);
}

void SO_Config_save_point::mpiReceive(int node) {
#ifdef USE_MPI
  MPI_Status status;
  int nelectrons;
  MPI_Recv(&nelectrons,1, MPI_INT, node, 0, MPI_Comm_grp,
           &status);
  electronpos.Resize(nelectrons);
  electronomega.Resize(nelectrons);
  for(int e=0; e< nelectrons; e++) {
    electronpos(e).Resize(3); 
    MPI_Recv(electronpos(e).v, 3, MPI_DOUBLE,
        node, 0, MPI_Comm_grp, &status);
    MPI_Recv(&electronomega(e), 1, MPI_DOUBLE,
        node, 0, MPI_Comm_grp, &status);

  }
#endif
}

void SO_Config_save_point::mpiSend(int node) {
#ifdef USE_MPI
  int nelectrons=electronpos.GetDim(0);
  MPI_Send(&nelectrons,1, MPI_INT, node, 0, MPI_Comm_grp);
  for(int e=0; e< nelectrons; e++) {
    MPI_Send(electronpos(e).v, 3, MPI_DOUBLE,
        node, 0, MPI_Comm_grp);
    MPI_Send(&electronomega(e), 1, MPI_DOUBLE,
        node, 0, MPI_Comm_grp);
  }
#endif
}




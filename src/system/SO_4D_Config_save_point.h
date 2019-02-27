#ifndef SO_4D_CONFIG_SAVE_POINT_H_INCLUDED
#define SO_4D_CONFIG_SAVE_POINT_H_INCLUDED
#include "Qmc_std.h"
#include "Sample_point.h"

class SO_4D_Config_save_point
{ 
 public:
  void write(ostream & os);
  void read(istream & is);
  void savePos(Sample_point * sample);
  void restorePos(Sample_point * smaple);
  void saveOmega(Sample_point * sample);
  void restoreOmega(Sample_point * sample);
  void getPos(int e, Array1 <doublevar> & r);
  void getOmega(int e, Array1 <doublevar> & omega);
  void mpiSend(int node);
  void mpiReceive(int node);
 private:
  Array1 < Array1 <doublevar > > electronpos;
  Array1 < Array1 <doublevar> > electronomega;
};

#endif

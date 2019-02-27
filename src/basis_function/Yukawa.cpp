

#include "Yukawa.h"
#include "qmc_io.h"

/*
 
*/
int Yukawa::read(
  vector <string> & words,
  unsigned int & pos)
{
  //cout << "gauss read " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  vector <string> alphatxt;



  pos=startpos;
  if(!readvalue(words, pos, expo1, "EXPO1") ) {
    error("Need EXPO1 in Yukawa");
  }
  pos=startpos;
  if(!readvalue(words, pos, expo2, "EXPO2") ) {
    error("Need EXPO2 in Yukawa");
  }

  return 0;
}

void Yukawa::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
//  parms.Resize(1);
  //parms(0)=log(gamma);
  //cout << "parms out " << parms(0) << endl;
}

void Yukawa::setVarParms(Array1 <doublevar> & parms) {

  //cout << "parms in " << parms(0) << endl;
}

int Yukawa::nfunc()
{
  return 1;
}

int Yukawa::showinfo(string & indent, ostream & os)
{
  os << indent << "Yukawa\n";
  os << indent << "EXPO1 " << expo1<<endl;
  os << indent << "EXPO2"  << expo2<<endl;
  return 1;
}

int Yukawa::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "YUKAWA\n";
  os << indent << "EXPO1 " << expo1 << endl;
  os << indent << "EXPO2 " << expo2 << endl;
  return 1;
}

void Yukawa::raw_input(ifstream & input)
{error("Raw input not supported by Yukawa");}

void Yukawa::calcVal(const Array1 <doublevar> & r,
                          Array1 <doublevar> & symvals,
                          const int startfill)
{
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= 1+startfill);
  
  symvals(startfill)=exp(-expo1*r(0))*(1-exp(-expo2*r(0)))*r(0);
  

}

void Yukawa::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  assert(r.GetDim(0) >=5);
  assert(symvals.GetDim(0) >= 1+startfill);
  assert(symvals.GetDim(1) >= 5);
  

  doublevar dudr,X, Y, dXdr, dYdr,
            d2udr2,drdx, drdy, drdz, d2rdx,d2rdy,d2rdz;
  drdx=r(2)/r(0); // dr/dx
  drdy=r(3)/r(0); // dr/dy
  drdz=r(4)/r(0); // dr/dz
  d2rdx=1/r(0)-r(2)*r(2)/(r(0)*r(1)); //d2r/dx2
  d2rdy=1/r(0)-r(3)*r(3)/(r(0)*r(1)); //d2r/dy2
  d2rdz=1/r(0)-r(4)*r(4)/(r(0)*r(1)); //d2r/dz2
  X=exp(-expo1*r(0)); //du/dr term1
  Y=-expo1-1/r(1)+((expo1+expo2)/r(0)+1/r(1))*exp(-expo2*r(0)); //du/dr term2
  dudr=X*Y;  // du/dr
  dXdr=-expo1*exp(-expo1*r(0));
  dYdr=2/(r(0)*r(1))+(-(expo1+expo2)/r(1)-2/(r(0)*r(1)))*exp(-expo2*r(0))-expo2*((expo1+expo2)/r(0)+1/r(1))*exp(-expo2*r(0));
  d2udr2=X*dYdr+dXdr*Y;
  symvals(startfill,0)=exp(-expo1*r(0))*(1-exp(-expo2*r(0)))*r(0);
 
  symvals(startfill,1)=dudr*drdx;
  symvals(startfill,2)=dudr*drdy;
  symvals(startfill,3)=dudr*drdz;
  symvals(startfill,4)=d2udr2*drdx*drdx+dudr*d2rdx
                      +d2udr2*drdy*drdy+dudr*d2rdy 
                      +d2udr2*drdz*drdz+dudr*d2rdz;
 
}

//------------------------------------------------------------------------

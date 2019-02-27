
#ifndef YUKAWA_H_INCLUDED
#define YUKAWA_H_INCLUDED

#include "Basis_function.h"

class Yukawa: public Basis_function
{
public:


  //-----------------------------------------------------


  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos

  );


  doublevar cutoff(int){ return 1; }

  int nfunc();
  virtual string label()
  {
    return centername;
  }

  int showinfo(string & indent, ostream & os);
  int writeinput(string &, ostream &);

  void raw_input(ifstream & input);

  void calcVal(const Array1 <doublevar> & r,
               Array1 <doublevar> & symvals,
               const int startfill=0);

  void calcLap(
    const Array1 <doublevar> & r,
    Array2 <doublevar> & symvals,
    const int startfill=0
  );
  virtual void calcHessian(const Array1 <doublevar> & r,
                           Array2 <doublevar> & symvals,
                           const int startfill=0){ }
  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);
  virtual int nparms() {
    return 1;
  }

private:

  string centername;
  doublevar expo1;
  doublevar expo2;
};

#endif // YUKAWA_H_INCLUDED
//--------------------------------------------------------------------------

/*
 
Copyright (C) 2007 Burkhard Militzer
 with some modifications by Lucas K. Wagner
 and extensions to Pfaffian algebra by Michal Bajdich

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
#ifndef MATRIXALGEBRA_H_INCLUDED
#define MATRIXALGEBRA_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"


/*!
  This is basically a replacement for doubles, but always with a logarithmic representation.
  I've tried to make it so it automatically does the right thing with standard 
  multiplication. Automatic downconversion to double is possible, but I've excluded it to 
  prevent accidental loss of precision.
  Note also that there may be performance problems with using these automatic conversions..
  */
struct log_real_value { 
  doublevar logval;
  int sign;
  doublevar val() { return sign*exp(logval); } 
  log_real_value() { logval=0; sign=1; }  //initializing to 1..
  //operator double() const { return sign*exp(logval); } 
  log_real_value(doublevar t) { 
    if(fabs(t) > 0)
      logval=log(fabs(t)); 
    else
      logval=-1e99;    //modified by MB 
    sign=t<0?-1:1; 
  }
};

inline log_real_value operator*(double t, log_real_value u) { 
  log_real_value v;
  v=u;
  if(fabs(t) > 0)
    v.logval+=log(fabs(t));
  else
    v.logval+=-1e99; //modified by MB 
      
  v.sign*=t<0?-1:1;
  return v;
}

inline log_real_value  operator*(log_real_value u,double t) { return t*u; }
inline log_real_value operator*(log_real_value t, log_real_value u) { 
  log_real_value v;
  v.logval=t.logval+u.logval;
  v.sign=t.sign*u.sign;
  return v;
}

//Try to safely sum a series of log_real_values
inline log_real_value sum(const Array1 <log_real_value> & vec) { 
  double s=0; 
  //for now, pivot on the first one
  //A more advanced way might be to check for underflows
  int n=vec.GetDim(0);
  int piv=0;
  for(int i=0; i< n; i++) { 
    s+=vec(i).sign*vec(piv).sign*exp(vec(i).logval-vec(piv).logval);
  }
  return s*vec(piv);
}

//--------------------------------------------------------------

doublevar dot(const Array1 <doublevar> & a,  const Array1 <doublevar> & b);
dcomplex  dot(const Array1 <doublevar> & a,  const Array1 <dcomplex> & b);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <doublevar> & b);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <dcomplex> & b);

//--------------------------------------------------------------

doublevar dot(const Array1 <doublevar> & a,  const Array1 <doublevar> & b, int & size);
dcomplex  dot(const Array1 <doublevar> & a,  const Array1 <dcomplex> & b, int & size);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <doublevar> & b, int & size);
dcomplex  dot(const Array1 <dcomplex> & a,  const Array1 <dcomplex> & b, int & size);


doublevar Determinant(const Array2 <doublevar> & a, const int n);
dcomplex Determinant(const Array2 <dcomplex> & a, const int n);

void MultiplyMatrices(const Array2 <doublevar> & a, const Array2 <doublevar> & b,
                      Array2 <doublevar> & c, int n);

void InvertMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                           const int lRow, const int n);
doublevar InverseUpdateRow(Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                           const int lRow, const int n);
doublevar InverseGetNewRatio(const Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
                             const int lCol, const int n);
doublevar InverseGetNewRatioRow(const Array2 <doublevar> & a1, const Array1 <doublevar> & newRow,
                             const int lRow, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                              const int lCol, const int n);
doublevar InverseUpdateColumn(Array2 <doublevar> & a1, const Array1 <doublevar> & newCol,
                              const int lCol, const int n);

void TransposeMatrix(Array2 <doublevar> & a, const int n);
log_real_value TransposeInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1, const int n);
doublevar TransposeInverseUpdateRow(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                    const int lCol, const int n);
doublevar TransposeInverseUpdateColumn(Array2 <doublevar> & a1, const Array2 <doublevar> & a,
                                       const int lRow, const int n);
doublevar Pfaffian_nopivot(const Array2 <doublevar> & a);

doublevar Pfaffian_partialpivot(const Array2 <doublevar> & a);
     
doublevar UpdateInversePfaffianMatrix(Array2 <doublevar> & a, Array1 <doublevar> & row,
                              Array1 <doublevar> & column, int n);
doublevar GetUpdatedPfaffianValue(Array2 <doublevar> & in, 
                                      Array1 <doublevar> & row, 
                                      int e);
doublevar PfaffianInverseMatrix(const Array2 <doublevar> & a, Array2 <doublevar> & a1);


inline doublevar cabs(dcomplex & a) {
  return sqrt(a.real()*a.real()+a.imag()*a.imag());
}

int ludcmp(Array2 <doublevar> & a, const int n, Array1 <int> & indx,
           doublevar & d);
void lubksb(Array2 <doublevar> & a, int n, Array1 <int> & indx,
            Array1 <doublevar> & b);

void EigenSystemSolverRealSymmetricMatrix(const Array2 < doublevar > & Ain, Array1 < doublevar> & evals, Array2 < doublevar> & evecs);

void GeneralizedEigenSystemSolverRealSymmetricMatrices(const Array2 < doublevar > & Ain, const Array2 < doublevar> & Bin, Array1 < doublevar> & evals, Array2 < doublevar> & evecs);

void sort_abs_values_descending(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, Array1 <int> & list);
void sort_abs_values_descending(const Array1 <dcomplex> & vals, Array1 <dcomplex> & newvals, Array1 <int> & list);
void sort_according_list(const Array1 <doublevar> & vals, Array1 <doublevar> & newvals, const Array1 <int> & list);
//------------------------------------------------------------------------
//IF Lapack was used
#ifdef USE_LAPACK
extern "C" { 
  void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
               double *A, int *LDAp, double *VLp, double *VUp,
               int *ILp, int *IUp, double *ABSTOLp, int *Mp,
               double *W, double *Z, int *LDZp, int *ISUPPZ,
               double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
               int *INFOp);
  double dlamch_(char *CMACHp);
  void dsygv_(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, 
              double *b, int *ldb, double *w, double *WORK, int *IWORK, int *info);
  void dgetrf_(int * m, int * n, double * A, int * lda, int * ipiv, int *info);
  void dgetrs_(char * trans, int *n, int * nrhs, double * A, int * lda, 
               int *ipiv, double * B, int * ldb, int * info);
};
#endif

//----------------------------------------------------------------------
//complex versions of determinant stuff

int ludcmp(Array2 <dcomplex > & a, const int n, 
           Array1 <int> & indx, doublevar & d);
void lubksb(Array2 < dcomplex > & a, int n, 
            Array1 <int> & indx, Array1 <dcomplex > & b);

dcomplex
TransposeInverseMatrix(const Array2 <dcomplex > & a, 
                       Array2 <dcomplex> & a1, 
                       const int n);
dcomplex
InverseUpdateColumn(Array2 <dcomplex > & a1, 
                    const Array1 <dcomplex > & newCol,
                    const int lCol, const int n);
void Jacobi(const Array2 < dcomplex > & Ain, Array1 <doublevar> & evals, Array2 < dcomplex > & evecs);

#endif // MATRIXALGEBRA_H_INCLUDED
//--------------------------------------------------------------------------

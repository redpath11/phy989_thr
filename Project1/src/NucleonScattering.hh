#ifndef NUCLEONSCATTERING_HH
#define NUCLEONSCATTERING_HH
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "lib.h"

//#define   ZERO       1.0E-10

class NucleonScattering
{
private:
  bool m_momentaSet;
  bool m_potentialSet;
  bool m_aSet;
  
  int m_meshPoints;
  double m_k0;
  double m_mass;
  
  double *m_weights;
  double *m_momenta;

  double **m_potential;
  double **m_A;
  double **m_R;

public:

  NucleonScattering(int,double);
  virtual ~NucleonScattering();

  double GetK0();
  void SetK0(double);

  void ListPoints();
  void PrintA();

  void GaussLegendreQuadrature(
    double x1,
    double x2,
    double *x,
    double *w,
    int n);
  void MapPoints();
  void FiniteSphere(double,double);
  void SetMatrixA();
  void PhaseShift();

};
#endif

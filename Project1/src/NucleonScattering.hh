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

  bool p_mapping;
  bool proximityWarning;
  double m_epsilon;
  
  int m_meshPoints;
  double m_k0;// [1/fm]
  double m_E0;// [1/fm]
  double m_V0,m_a;// [1/fm],[fm]
  double m_mass;// [1/fm]
  
  double *m_weights;
  double *m_momenta;

  double **m_potential;
  double **m_A;
  double **m_R;

public:

  NucleonScattering(int,double);
  ~NucleonScattering();

  double GetE0();
  void SetE0(double);
  double GetK0();
  void SetK0(double);

  double GetV0();
  void SetV0(double);

  double GetA();
  void SetA(double);

  void SetProxWarning(){
    proximityWarning=true;
  }
  bool CheckProxWarning(){
    return proximityWarning;
  }

  double GetEpsilon(){
    return m_epsilon;
  }

  void ListPoints();
  void PrintV();
  void PrintA();
  void PrintIdentityHopefully();

  void GaussLegendreQuadrature(
    double x1,
    double x2,
    double *x,
    double *w,
    int n);
  void MapPoints();
  void FiniteSphere(double,double);
  void PotentialNP();
  void SetMatrixA();
//  void PhaseShift();
  double PhaseShift();

  double AnalyticPhaseShift();

};
#endif

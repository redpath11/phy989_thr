//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function */
int main()
{
/**/
  NucleonScattering n(100,5.);
  n.MapPoints();
//  n.SetK0(1.2);
//  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
  n.FiniteSphere(0.5,1.);
//  n.PotentialNP();
  n.SetMatrixA();
//  n.PrintA();
  n.PhaseShift();
  n.AnalyticPhaseShift();
//  n.~NucleonScattering();

/*
  NucleonScattering n(10,2.0);
  n.SetV0(0.5);
  n.SetA(1.);
  n.AnalyticPhaseShift();
*/
}

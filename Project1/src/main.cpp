//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function */
int main()
{
/**/
  NucleonScattering n(1000,1.);
  n.MapPoints();
  n.SetK0(1.2);
//  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
  n.FiniteSphere(45.,0.5);
  n.SetMatrixA();
//  n.PrintA();
  n.PhaseShift();
  n.AnalyticPhaseShift();
//  n.~NucleonScattering();

/*
  NucleonScattering n(10,1.0);
  n.SetV0(45.);
  n.SetA(0.5);
  n.AnalyticPhaseShift();
*/
}

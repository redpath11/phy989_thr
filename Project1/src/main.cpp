//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function */
int main()
{
  NucleonScattering n(10,1.);
  n.MapPoints();
  n.SetK0(1.2);
  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
  n.FiniteSphere(1.,2.);
  n.SetMatrixA();
//  n.PrintA();
  n.PhaseShift();
//  n.~NucleonScattering();
}

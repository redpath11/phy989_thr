//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- single energy version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab = atof(argv[2]);
  NucleonScattering n(nPts,eLab);
  n.MapPoints();
//  n.SetK0(eLab);
//  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
  n.FiniteSphere(-0.5,1.);
//  n.PotentialNP();
  n.SetMatrixA();
  n.PrintV();
  n.PhaseShift();
//  n.PrintA();
  n.AnalyticPhaseShift();
//  n.~NucleonScattering();

/*
  NucleonScattering n(10,2.0);
  n.SetV0(0.5);
  n.SetA(1.);
  n.AnalyticPhaseShift();
*/
}

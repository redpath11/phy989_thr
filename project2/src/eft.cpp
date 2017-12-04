//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- single energy version */
int main(int argc, char *argv[])
{
  if(argc!=7){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab = atof(argv[2]);
  double eftC0 = atof(argv[3]);
  double eftC2 = atof(argv[4]);
  double eftC4 = atof(argv[5]);
  double eftC4p= atof(argv[6]);

  NucleonScattering n(nPts,eLab);
  n.MapPoints();
//  n.FiniteSphere(-0.5,1.);
  n.PotentialNP();
  n.SetMatrixA();
//  n.PrintV();
  cout << "NP delta: " << n.PhaseShift() << endl;

  NucleonScattering p(nPts,eLab);
  p.MapPoints();
//  p.eftLO(eftC0);
//  p.eftNLO(eftC0,eftC2);
  p.eftNNLO(eftC0,eftC2,eftC4,eftC4p);
  p.SetMatrixA();
  cout << "EFT delta: " << p.PhaseShift() << endl;

//  n.PrintA();
//  n.AnalyticPhaseShift();
//  n.~NucleonScattering();

/*
  NucleonScattering n(10,2.0);
  n.SetV0(0.5);
  n.SetA(1.);
  n.AnalyticPhaseShift();
*/
}

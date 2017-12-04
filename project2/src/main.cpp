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
  ofstream output("../output/eftLO.txt");
//  ofstream output("../output/eftNLO.txt");
//  ofstream output("../output/eftNNLO.txt");
  eLab = 0.;

  for(int i=0;i<25;i++){
    eLab += 0.5;

  // emperical NP potential
  NucleonScattering n(nPts,eLab);
  n.MapPoints();
//  n.FiniteSphere(-0.5,1.);
  n.PotentialNP();
  n.SetMatrixA();
//  n.PrintV();
  output << setw(12) << eLab
  	<< setw(12) << n.PhaseShift();

  // EFT potential
  NucleonScattering p(nPts,eLab);
  p.MapPoints();
  p.eftLO(-0.985);
//  p.eftNLO(-0.985,0.75);
//  p.eftNNLO(-0.985,0.75,0.2,0.2);
  p.SetMatrixA();

  output << setw(12) << p.PhaseShift() << endl;

  }// close loop over energies

  output.close();

}

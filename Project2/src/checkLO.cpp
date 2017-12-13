//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- loop over energies version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]> <CO>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab = atof(argv[2]);
  double eftC0 = atof(argv[3]);

  ofstream output("../output/LO.txt");

//  for(int i=0;i<50;i++){
//  for(int i=0;i<20;i++){
  for(int i=0;i<15;i++){
    NucleonScattering n(nPts,eLab);
    n.MapPoints();
    n.PotentialNP();
    n.SetMatrixA();
    output << setw(12) << eLab
         << setw(12) << n.PhaseShift();

    // EFT potential
    NucleonScattering p(nPts,eLab);
    p.MapPoints();
    p.eftLO(eftC0);
//    p.OnePionLO(eftC0);
    p.SetMatrixA();
    output << setw(12) <<  p.PhaseShift() << endl;

//    eLab += 0.1;
    eLab *= 2.;

  }

  output.close();
}

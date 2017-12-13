//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- loop over energies version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]> <C0> <C2>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab = atof(argv[2]);
  double eftC0 = atof(argv[3]);
  double eftC2 = atof(argv[4]);

  ofstream output("../output/NLO.txt");

//  for(int i=0;i<50;i++){
  for(int i=0;i<20;i++){
    NucleonScattering n(nPts,eLab);
    n.MapPoints();
    n.PotentialNP();
    n.SetMatrixA();
    output << setw(12) << eLab
         << setw(12) << n.PhaseShift();

    // EFT potential
    NucleonScattering p(nPts,eLab);
    p.MapPoints();
//    p.eftNLO(eftC0,eftC2);
    p.OnePionNLO(eftC0,eftC2);
    p.SetMatrixA();
    output << setw(12) <<  p.PhaseShift() << endl;

//    eLab += 0.1;
    eLab *= 2.0;

  }

  output.close();
}

//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- loop over energies version */
int main(int argc, char *argv[])
{
  if(argc!=6){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]> <C0> <C2> <C4> <C4p>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab = atof(argv[2]);
  double eftC0 = atof(argv[3]);
  double eftC2 = atof(argv[4]);
  double eftC4 = atof(argv[5]);
  double eftC4p = atof(argv[6]);

  ofstream output("../output/NNLO.txt");

//  for(int i=0;i<50;i++){
  for(int i=0;i<20;i++){
    NucleonScattering n(nPts,eLab);
    n.MapPoints();
    n.PotentialNP();
    n.SetMatrixA();
//    double np = n.PhaseShift();
    output << setw(12) << eLab
         << setw(12) << n.PhaseShift();
//         << setw(12) << np;

    // EFT potential
    NucleonScattering p(nPts,eLab);
    p.MapPoints();
//    p.eftNNLO(eftC0,eftC2,eftC4,eftC4);
    p.OnePionNNLO(eftC0,eftC2,eftC4,eftC4p);
    p.SetMatrixA();
//    double eft = p.PhaseShift();
    output << setw(12) <<  p.PhaseShift() << endl;
//    output << setw(12) <<  eft << endl;

//    eLab += 0.1;
    eLab *= 2.0;

  }

  output.close();
}

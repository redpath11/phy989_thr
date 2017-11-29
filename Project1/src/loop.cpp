//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- loop over energies version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy [MeV]>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  double eLab[11] = {1.,5.,10.,25.,50.,100.,150.,200.,250.,300.,350.};
  ofstream output("../output/l0ps.txt");

  for(int i=0;i<11;i++){
    NucleonScattering n(nPts,eLab[i]);
    n.MapPoints();
//  n.SetK0(1.2);
//  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
//    n.FiniteSphere(0.5,1.);
    n.PotentialNP();
    n.SetMatrixA();
//  n.PrintA();
    output << setw(12) << eLab[i]
         << setw(12) << n.PhaseShift() << endl;

  }

  output.close();
}

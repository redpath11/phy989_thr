//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- loop over energies version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <# mesh points> <lab energy>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
//  double eLab[11] = {1.,5.,10.,25.,50.,100.,150.,200.,250.,300.,350.};
//  double eLab[11] = {1.,2.,5.,8.,10.,15.,17.,25.,30.,40.,50.};
//  ofstream output("../output/l0ps.txt");
  ofstream output("../output/FiniteSpherel0ps.txt");
  double eLab = 0.;

  for(int i=0;i<25;i++){
    eLab += 0.5;
    NucleonScattering n(nPts,eLab);
    n.MapPoints();
//  n.SetK0(1.2);
//  n.ListPoints();
//  cout << endl << "pi/4 = " << M_PI_4 << endl;
    n.FiniteSphere(-0.5,1.);
//    n.PotentialNP();
    n.SetMatrixA();
//  n.PrintA();
    output << setw(12) << n.GetE0()
//    	 << setw(12) << n.GetK0()
         << setw(12) << n.PhaseShift()
	 << setw(12) << n.AnalyticPhaseShift() << endl;

  }

  output.close();
}

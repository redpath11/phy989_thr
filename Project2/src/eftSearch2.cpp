//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- C0,C2 search version */
int main(int argc, char *argv[])
{
  if(argc!=5){
    cout << "Usage: ./ns.exe <# mesh points> <startC0> <incC0> <startC2> <incC2>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  const double searchC0 = atof(argv[2]);
  const double incC0 = atof(argv[3]);
  const double searchC2 = atof(argv[4]);
  const double incC2 = atof(argv[5]);

  double eftC0 = searchC0;
  double eftC2 = searchC2;
  ofstream output("../output/eftNLOfit.txt");


for(int i=0;i<50;i++){// loop over values of eftC0
  for(int j=0;j<50;j++){// loop over values of eftC2

    double eLab = 0.001;
    double weight = 1.;
    double chi2 = 0.;

//    for(int k=0;k<5;k++){// loop over energies
      // emperical NP potential
      NucleonScattering n(nPts,eLab);
      n.MapPoints();
      n.PotentialNP();
      n.SetMatrixA();
      double np = n.PhaseShift();
      
      // EFT potential
      NucleonScattering p(nPts,eLab);
      p.MapPoints();
      p.eftNLO(eftC0,eftC2);
//      p.OnePionNLO(eftC0,eftC2);
      p.SetMatrixA();
      double ef = p.PhaseShift();

      chi2 += pow((np - ef),2.)/weight;// weight low energies
      eLab += 0.001;// increment eLab
      weight *= 2.;

//    }
//      chi2 *= 0.2;// all energies weighted evenly
      output << setw(12) << eftC0
		<< setw(12) << eftC2
		<< setw(12) << chi2 << endl;
  eftC2 += incC2;// increment eftC2
  }// close loop over eftC2
  eftC0 += incC0;// increment eftC0
  eftC2 = searchC2;
}// close loop over eftC0

  return 0;

}

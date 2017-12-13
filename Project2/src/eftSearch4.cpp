//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- C0,C2 search version */
int main(int argc, char *argv[])
{
  if(argc!=9){
    cout << "Usage: ./ns.exe <# mesh points> <startC0> <incC0> <startC2> <incC2> <startC4> <incC4>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  const double startC0 = atof(argv[2]);
  const double incC0 = atof(argv[3]);
  const double startC2 = atof(argv[4]);
  const double incC2 = atof(argv[5]);
  const double startC4 = atof(argv[6]);
  const double incC4 = atof(argv[7]);
  const double startC4p = 0.;// atof(argv[8]);
  const double incC4p = 0.;//atof(argv[9]);

  double eftC0 = startC0;
  double eftC2 = startC2;
  double eftC4 = startC4;
  double eftC4p = startC4p;
  ofstream output("../output/eftNNLOfit.txt");


for(int i=0;i<5;i++){// loop over values of eftC0
  for(int j=0;j<5;j++){// loop over values of eftC2
    for(int k=0;k<5;k++){// loop over values of eftC4
//      for(int l=0;l<5;l++){// loop over values of eftC4p

    double eLab = 0.001;
    double chi2 = 0.;
    double weight = 1.;

    for(int energy=0;energy<5;energy++){// loop over energies
      // emperical NP potential
      NucleonScattering n(nPts,eLab);
      n.MapPoints();
      n.PotentialNP();
      n.SetMatrixA();
      double np = n.PhaseShift();
      
      // EFT potential
      NucleonScattering p(nPts,eLab);
      p.MapPoints();
//      p.eftNNLO(eftC0,eftC2,eftC4,eftC4p);
      p.OnePionNNLO(eftC0,eftC2,eftC4,eftC4p);
      p.SetMatrixA();
      double ef = p.PhaseShift();

      chi2 += pow((np - ef),2.)/weight;// weight low energies
      eLab += 0.001;// increment eLab
      weight *= 2.;

    }
//      chi2 *= 0.2;// all energies weighted evenly
      output << setw(12) << eftC0
		<< setw(12) << eftC2
		<< setw(12) << eftC4
		<< setw(12) << eftC4p
		<< setw(12) << chi2 << endl;

      eftC4p *= incC4p;
//      }// close loop over eftC4p
    eftC4 *= incC4;
    eftC4p = startC4p;
    }// close loop over eftC4
  eftC2 *= incC2;// increment eftC2
  eftC4 = startC4;
  eftC4p = startC4p;
  }// close loop over eftC2
  eftC0 *= incC0;// increment eftC0
  eftC2 = startC2;
  eftC4 = startC4;
  eftC4p = startC4p;
}// close loop over eftC0

  return 0;

}

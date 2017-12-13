//#include "lib.h"
#include "NucleonScattering.hh"

/* The Main Function -- single energy version */
int main(int argc, char *argv[])
{
  if(argc!=3){
    cout << "Usage: ./ns.exe <starting C0> <C0 increment>" << endl;
  }
/**/
  int nPts = atoi(argv[1]);
  const double searchC0 = atof(argv[2]);
  const double incC0 = atof(argv[3]);
  ofstream output("../output/eftLOfit.txt");

  double eftC0 = searchC0;


for(int h=0;h<100;h++){// loop over values of eftC0
  double eLab = 0.001;
  double weight = 1.;
  double chi2 = 0.;

  for(int i=0;i<5;i++){// loop over low energies

  // emperical NP potential
  NucleonScattering n(nPts,eLab);
  n.MapPoints();
  n.PotentialNP();
  n.SetMatrixA();
  double np = n.PhaseShift();

  // EFT potential
  NucleonScattering p(nPts,eLab);
  p.MapPoints();
  p.eftLO(eftC0);
//  p.OnePionLO(eftC0);
  p.SetMatrixA();
  double ef = p.PhaseShift();

  chi2 += pow((np - ef),2.)/weight;// weight low energies more
  eLab += 0.001;// increment eLab
  weight *= 2.;// adjust weight

  }// CLOSE loop over low energies
  
//  chi2 *= 0.2;// all energies weighted evenly
  output << setw(12) << eftC0
  	<< setw(12) << chi2 << endl;
//  printf("eftC0 = %g   chiSq = %g\n",eftC0,chi2);
  eftC0 += incC0;// increment eftC0
}

  return 0;

}

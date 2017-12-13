//#include "lib.h"
#include "NucleonScattering.hh"

double calculateX2(int np,double C0,double C2,double C4,double C4p);

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
  const double startC4p = atof(argv[8]);
  const double incC4p = atof(argv[9]);

  double eftC0 = startC0;
  double eftC2 = startC2;
  double eftC4 = startC4;
  double eftC4p = startC4p;
  ofstream output("../output/eftNNLOfit.txt");
  int nIterations=0;
  double startX2 = 1.;

  while(nIterations < 1e3 && startX2>1e-6){
  // calculate chi2 at the start point
  startX2 = calculateX2(nPts,eftC0,eftC2,eftC4,eftC4p);
  double next[8]  = {0.,0.,0.,0.,0.,0.,0.,0.};
  // calculate chi2 at C0 +/- C0inc
  next[0] = calculateX2(nPts,eftC0-incC0,eftC2,eftC4,eftC4p);
  next[1] = calculateX2(nPts,eftC0+incC0,eftC2,eftC4,eftC4p);
  // calculate chi2 at C2 +/- C2inc
  next[2] = calculateX2(nPts,eftC0,eftC2-incC2,eftC4,eftC4p);
  next[3] = calculateX2(nPts,eftC0,eftC2+incC2,eftC4,eftC4p);
  // calculate chi2 at C4 +/- C4inc
  next[4] = calculateX2(nPts,eftC0,eftC2,eftC4-incC4,eftC4p);
  next[5] = calculateX2(nPts,eftC0,eftC2,eftC4+incC4,eftC4p);
  // calculate chi2 at c4p +/- C4pinc
  next[6] = calculateX2(nPts,eftC0,eftC2,eftC4,eftC4p-incC4p);
  next[7] = calculateX2(nPts,eftC0,eftC2,eftC4,eftC4p+incC4p);

  output << setw(12) << eftC0
  	<< setw(12) << eftC2
	<< setw(12) << eftC4
	<< setw(12) << eftC4p
	<< setw(12) << startX2 << endl;

  double min = next[0];
  int index = 0;

  // find direction of steepest X2 change
  for(int i=0;i<8;i++){
    if((next[i] < startX2) && (next[i] < min)){
      min = next[i];
      index = i;
    }
  }

  if(index==0){eftC0 -= incC0;}
  else if(index==1){eftC0 += incC0;}
  else if(index==2){eftC2 -= incC2;}
  else if(index==3){eftC2 += incC2;}
  else if(index==4){eftC4 -= incC4;}
  else if(index==5){eftC4 += incC4;}
  else if(index==6){eftC4p -= incC4p;}
  else if(index==7){eftC4p += incC4p;}

  printf("New Coordintaes %g,%g,%g,%g\n",eftC0,eftC2,eftC4,eftC4p);
  nIterations ++;

  }// close while

  printf("Minimum Coordintaes %g,%g,%g,%g\n",eftC0,eftC2,eftC4,eftC4p);

  output.close();

  return 0;

}

double calculateX2(int np,double C0,double C2,double C4,double C4p){
  double X2=0.;
  double eLab = 0.001;
  double weight = 1.;

  for(int energy=0;energy<5;energy++){// loop over energies
    // emperical NP potential
    NucleonScattering n(np,eLab);
    n.MapPoints();
    n.PotentialNP();
    n.SetMatrixA();
    double np = n.PhaseShift();

    // EFT potential
    NucleonScattering p(np,eLab);
    p.MapPoints();
//      p.eftNNLO(C0,C2,C4,C4p);
    p.OnePionNNLO(C0,C2,C4,C4p);
    p.SetMatrixA();
    double ef = p.PhaseShift();

    X2 += pow((np - ef),2.)/weight;// weight low energies
    eLab += 0.001;// increment eLab
    weight *= 2.;

  }

  return X2;

}

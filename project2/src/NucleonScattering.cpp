#include "NucleonScattering.hh"

NucleonScattering::NucleonScattering(int n,double Elab=1.){
  m_momentaSet=false; m_potentialSet=false; m_aSet=false;
  m_meshPoints=n;

  // Set the lengths for the k and w arrays
  m_weights = new double [m_meshPoints+1];
  m_momenta = new double [m_meshPoints+1];

  // Set the nucleon mass
  m_mass = 939./(197.);// [1/fm]
//  m_mass = 939.;// [MeV]
//  m_mass = 1.0;// [Mn]

  // Convert the lab energy to [1/fm] units and compute k0
  m_E0 = Elab/197.;// [1/fm]
  m_k0 = sqrt(Elab / 83.);// [1/fm]


  // Set flags to print stuff out
  p_mapping = false;

  m_epsilon = 0.00001;
  proximityWarning=false;

  m_Lambda = 0.343;// ~ 1/2 pion mass
}

NucleonScattering::~NucleonScattering(){

  delete m_weights; delete m_momenta;
//  free_matrix((void **) m_potential);
};


// Functions to set and get values for k0 and E0
double NucleonScattering::GetK0(){ return m_k0; }
void NucleonScattering::SetK0(double k0){
  // set observable point
  m_k0 = k0;
  m_momenta[m_meshPoints] = m_k0;
  m_E0 = m_k0*m_k0*83.;// Elab
}

double NucleonScattering::GetV0(){ return m_V0; }
void NucleonScattering::SetV0(double v){ m_V0 = v; }

double NucleonScattering::GetA(){ return m_a; }
void NucleonScattering::SetA(double a){ m_a = a; }

double NucleonScattering::GetE0(){ return m_E0; }
void NucleonScattering::SetE0(double e0){
  m_E0 = e0/197.;// [1/fm]
  m_k0 = sqrt(m_E0/83.);
  m_momenta[m_meshPoints] = m_k0;
  cout << "E0 = " << e0   << " MeV" << endl;
  cout << "   = " << m_E0 << " 1/fm" << endl;
  cout << "k0 = " << m_k0 << " 1/fm" << endl;
}



       /* stolen from the libraries provided
       **
       ** The function
       **              GaussLegendreQuadrature()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

//void NucleonScattering::GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n)
void NucleonScattering::GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
           ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
         p2 =0.0;

           /*
           ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

         for(j = 1; j <= n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
         }

           /*
           ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

         pp = n * (z * p1 - p2)/(z * z - 1.0);
         z1 = z;
         z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
          ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function GaussLegendreQuadrature()


/* This function sets up the mesh points and weights
 * and then maps them from the interval [-1,1] to the
 * interval [0,inf). Note that the units of momentum
 * are [fm^-1].
*/
void NucleonScattering::MapPoints(){
  // Call function to fill the weights and k points arrays

  double *weights = new double [m_meshPoints];
  double *momenta = new double [m_meshPoints];

  GaussLegendreQuadrature(-1.,1.,momenta,weights,m_meshPoints);

  // Now map w and k from [-1,1] to [0,inf)
  for(int i=0;i<m_meshPoints;i++){
  if(p_mapping){ printf("x = %g   w = %g",momenta[i],weights[i]);}
    m_momenta[i] = tan(M_PI_4*(1.+momenta[i]));
//    m_momenta[i] = 197.*tan(M_PI_4*(1.+momenta[i]));
    double costerm = cos(M_PI_4*(1.+momenta[i]));
    m_weights[i] = M_PI_4 * weights[i] / (costerm*costerm);
//    m_weights[i] = 197.*M_PI_4 * weights[i] / (costerm*costerm);
  if(p_mapping){  printf("    k = %g   w = %g\n",m_momenta[i],m_weights[i]);}
  }

  m_momenta[m_meshPoints] = m_k0;
  m_weights[m_meshPoints] = 0.;

  m_momentaSet=true;
}


void NucleonScattering::ListPoints(){
  for(int i=0;i<m_meshPoints+1;i++){
    printf("k[%i] = %g\n",i,m_momenta[i]);
  }
}


/* This function sets up the potential matrix for a
 * finite sphere of radius a and strength V
*/
void NucleonScattering::FiniteSphere(double V,double a){
m_V0 = V;
m_a  = a;// [fm]
if(m_momentaSet){
  // allocate memory for potential matrix
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_potential[i][j] = 0.;
    }
  }

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<i;j++){
      if(m_momenta[i]!=m_momenta[j]){
      double num1 = m_momenta[j] * sin(m_momenta[i] * m_a) * cos(m_momenta[j] * m_a);
      double num2 = m_momenta[i] * cos(m_momenta[i] * m_a) * sin(m_momenta[j] * m_a);
      double denom= pow(m_momenta[i],3.)*m_momenta[j] - m_momenta[i]*pow(m_momenta[j],3.);
      m_potential[i][j] = m_V0 * (num1 - num2) / denom;
      }
      else{
        m_potential[i][j] = 0.5*m_a - (sin(2.*m_momenta[i]*m_a)/(4.*m_momenta[i]));
	m_potential[i][j] *= (m_V0 / pow(m_momenta[i],2.));
      }
      m_potential[j][i] = m_potential[i][j];

    }// Close loop over j
    m_potential[i][i] = 0.5*m_a - (sin(2.*m_momenta[i]*m_a)/(4.*m_momenta[i]));
    m_potential[i][i] *= (m_V0 / pow(m_momenta[i],2.));
  }

  m_potentialSet=true;
}
else if(m_potentialSet){cout << "Potential already set." << endl;return;}
else{cout << "Run MapPoints() first." << endl;return;}
}


/* Function to set up the NP potential as a sum of 3 Yukawa
 * potentials with ranges 1, 4, 7 times the inverse pion mass.
 * The well depths Va, Vb, Vc are given in units of MeV, we
 * convert them to unitless factors
*/ 
void NucleonScattering::PotentialNP(){
  const double u = 0.7;// [1/fm]
//  const double u = 0.7*197.;// [MeV]
  double Va = -10.463/(197.*u);// []
//  double Va = -10.463/(u);// []
  double  a =  1.*u;
  a*=a;
  double Vb = -1650.6/(197.*u);// []
//  double Vb = -1650.6/(u);// []
  double  b = 4.*u;
  b*=b;
  double Vc =  6484.3/(197.*u);// []
//  double Vc =  6484.3/(u);// []
  double  c = 7.*u;
  c*=c;

if(m_momentaSet){
  // allocate memory for potential matrix
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_potential[i][j] = 0.;
    }
  }

  cout << m_potential << endl;
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<=i;j++){
      double potA = 0.;
      potA       += log(pow(m_momenta[i]+m_momenta[j],2.) + a);
      potA       -= log(pow(m_momenta[i]-m_momenta[j],2.) + a);
      potA       *= Va / (4. * m_momenta[i]*m_momenta[j]);// [fm^2]
/**/
      double potB = log(pow(m_momenta[i]+m_momenta[j],2.) + b);
      potB       -= log(pow(m_momenta[i]-m_momenta[j],2.) + b);
      potB       *= Vb / (4. * m_momenta[i]*m_momenta[j]);// [fm^2]

      double potC = log(pow(m_momenta[i]+m_momenta[j],2.) + c);
      potC       -= log(pow(m_momenta[i]-m_momenta[j],2.) + c);
      potC       *= Vc / (4. * m_momenta[i]*m_momenta[j]);// [fm^2]

      m_potential[i][j] = potA + potB + potC;
      m_potential[j][i] = m_potential[i][j];
//      printf("  V[%i][%i] = %g",i,j,m_potential[i][j]);
    }
//      cout << endl;
  }
  m_potentialSet=true;
}
else if(m_potentialSet){cout << "Potential already set." << endl;return;}
else{cout << "Run MapPoints() first." << endl;return;}
}


/* Print the potential matrix to a file */
void NucleonScattering::PrintV(){
if(m_potentialSet){
ofstream vOut("../output/checkPotential.txt");
  for(int i=0;i<=m_meshPoints;i++){
    for(int j=0;j<=m_meshPoints;j++){
//      printf("  V[%i][%i] = %g",i,j,m_potential[i][j]);
      vOut << setw(15) << i << "  " 
           << setw(15) << j << "  "
	   << setw(15) << m_potential[i][j] << endl;
    }
  }
vOut.close();
}
else {cout << "[V] not set." << endl;return;}
}


void NucleonScattering::eftLO(double C0){
if(m_momentaSet){
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_potential[i][j] = exp(-1.*pow(m_momenta[i],4.)/pow(m_Lambda,4.));
      m_potential[i][j] *= C0;
      m_potential[i][j] *= exp(-1.*pow(m_momenta[j],4.)/pow(m_Lambda,4.));
    }
  }
  m_potentialSet=true;
}
}


void NucleonScattering::eftNLO(double C0,double C2){
if(m_momentaSet){
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_potential[i][j] = exp(-1.*pow(m_momenta[i],4.)/pow(m_Lambda,4.));
      m_potential[i][j] *= (C0 + C2*(pow(m_momenta[i],2.) + pow(m_momenta[j],2.)));
      m_potential[i][j] *= exp(-1.*pow(m_momenta[j],4.)/pow(m_Lambda,4.));
    }
  }
  m_potentialSet=true;
}
}


void NucleonScattering::eftNNLO(double C0,double C2,double C4,double C4p){
if(m_momentaSet){
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      double Vij = C0;
      Vij += C2*(pow(m_momenta[i],2.) + pow(m_momenta[j],2.));
      Vij += C4*(pow(m_momenta[i],4.) + pow(m_momenta[j],4.));
      Vij += C4p*pow(m_momenta[i],2.)*pow(m_momenta[j],2.);
      m_potential[i][j] = exp(-1.*pow(m_momenta[i],4.)/pow(m_Lambda,4.));
      m_potential[i][j] *= Vij;
      m_potential[i][j] *= exp(-1.*pow(m_momenta[j],4.)/pow(m_Lambda,4.));
    }
  }
  m_potentialSet=true;
}
}


/* Compute the matrix elements A(i,j). There are four cases to
 * consider :
 * (1) A(i,j) = DiracDelta(i,j) - V(i,j) * u(j); i,j from 1..N
 * (2) A(N+1,j) = -V(i,j) * u(j); j from 1..N (last row)
 * (3) A(i,N+1) = -V(i,N+1) * u(N+1); i from 1..N (last column)
 * (4) A(N+1,N+1) = 1 - V(N+1,N+1) * u(N+1) (bottom right corner)
*/
void NucleonScattering::SetMatrixA(){
if(m_potentialSet){
  cout << "Setting matrix A." << endl;
  const double TwoOverPi = M_2_PI;
  // allocate memory for [A]
  m_A = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  cout << m_A << endl;
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_A[i][j] = 0.;
    }
  }

  // calculate u(N+1)
  double u0 = 0.;
  for(int i=0;i<m_meshPoints;i++){
    u0 += m_weights[i]/(m_k0*m_k0 - pow(m_momenta[i],2.));
  }
  u0 *= (-1.)*TwoOverPi * m_mass * m_k0 * m_k0;
  // now we have u(N+1)

  // compute for cases (1) and (2)
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints;j++){
      // Calculate uj terms
      double num   = m_mass* m_weights[j] * m_momenta[j]*m_momenta[j];
      double denom = (m_k0*m_k0 - pow(m_momenta[j],2.));
      double uj = TwoOverPi * num / denom;
      // Calculate elements of [A]
      m_A[i][j] = (-1.)*m_potential[i][j]*uj;
    }
    m_A[i][i] += 1.;
  }

  // case (3)
  for(int i=0;i<m_meshPoints;i++){
    m_A[i][m_meshPoints] = (-1.)*m_potential[i][m_meshPoints]*u0;
/*    
    double epsilon = abs(pow(m_momenta[j],2.) - pow(m_k0,2.));
    if(epsilon<m_epsilon){
      printf("WARNING: k0 close to a mesh point: %g\n",epsilon);
      SetProxWarning();
    }
*/
    
  }

  // case (4) - the last element of [A]
  m_A[m_meshPoints][m_meshPoints] = 1. - m_potential[m_meshPoints][m_meshPoints]*u0;



//  m_A[m_meshPoints][m_meshPoints] = 1. - m_potential[m_meshPoints][m_meshPoints]*u0;
  m_aSet=true;
}
}


void NucleonScattering::PrintA(){
if(m_aSet){
ofstream vOut("../output/checkA.txt");
  for(int i=0;i<=m_meshPoints;i++){
    for(int j=0;j<=m_meshPoints;j++){
      vOut << setw(15) << i << "  " 
           << setw(15) << j << "  "
	   << setw(15) << m_A[i][j] << endl;
    }
  }
vOut.close();
}
else {cout << "[A] not set." << endl;return;}
}


void NucleonScattering::PrintIdentityHopefully(){
if(m_aSet){
  ofstream iOut("../output/checkIdentity.txt");
  double **aInverse = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      aInverse[i][j] = m_A[i][j];
    }
  }
  inverse(aInverse,m_meshPoints+1);
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      double cell = 0.;
      for(int k=0;k<m_meshPoints+1;k++){
        cell += m_A[i][k]*aInverse[k][j];
      }
      iOut << setw(18) << i
           << setw(18) << j
	   << setw(18) << cell << endl;;
    }
  }
  iOut.close();
}
else{cout << "Set A first." << endl;}
}


/* Calculate the phase shift from R(N+1,N+1) */

//void NucleonScattering::PhaseShift(){
double NucleonScattering::PhaseShift(){
if(m_aSet){
  // invert [A]
  inverse(m_A,m_meshPoints+1);
  // allocate memory for [R] and initialize
  m_R = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));

  // perform matrix multiplication
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_R[i][j] = 0.;
      for(int k=0;k<m_meshPoints+1;k++){
        m_R[i][j] += m_A[i][k] * m_potential[k][j];
      }
    }
  }

//  printf("R(k0,k0) = %g\n",m_R[m_meshPoints][m_meshPoints]);
  double delta = (-1.0)*m_R[m_meshPoints][m_meshPoints]*m_mass*m_k0;
//  printf("tan(delta) = %g\n",delta);
  delta = atan(delta);
//  printf("delta0 = %g\n",delta);
  /*
  if(delta<0){
    delta += M_PI;
    printf("delta1 = %g\n",delta);
  }
  */
//  delta = atan(tan(delta));
//  printf("atan(tan(delta)) = %g\n",atan(tan(delta)));
//  printf("delta1 = %g degrees\n",delta*180./M_PI);
  return delta;
//  return atan(tan(delta));
}
else{cout << "Run SetMatrixA() first." << endl;return 0.;}
}


double NucleonScattering::AnalyticPhaseShift(){
// for equal mass scattering, the reduced mass is
// 0.5 * the particle mass
  double v0 = -1. * m_V0;
//  double E = (m_k0*m_k0) / (m_mass);
  double delta = sqrt(m_E0/(m_E0+v0));
  delta *= tan(m_a*sqrt(m_mass*(m_E0+v0)));
  delta = atan(delta);
  if(delta<0){delta += M_PI;}
  delta -= m_a*sqrt(m_mass*m_E0);

//  printf("mu = %g\na = %g\nV0 = %g\n",reducedMass,m_a,m_V0);
//  printf("k = %g\n",m_k0);
  printf("E = %g\nk0 = %g\n",m_E0,m_k0);
//  printf("R(k0,k0)= %g\n",-1.*tan(delta)/(m_mass*m_k0));
  printf("AnalyticDelta = %g\n",delta);

//  return atan(tan(delta));
  return delta;
}

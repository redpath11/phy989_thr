#include "NucleonScattering.hh"

NucleonScattering::NucleonScattering(int n,double k0=1.){
  m_momentaSet=false; m_potentialSet=false; m_aSet=false;
  m_meshPoints=n;
  m_k0 = k0;
  m_weights = new double [m_meshPoints+1];
  m_momenta = new double [m_meshPoints+1];
//  m_mass = 940.;// [MeV/c^2]
  m_mass = 1.;// [MeV/c^2]

}

NucleonScattering::~NucleonScattering(){

  delete m_weights; delete m_momenta;
//  free_matrix((void **) m_potential);
};

double NucleonScattering::GetK0(){ return m_k0; }
void NucleonScattering::SetK0(double k0){
  // set observable point
  m_k0 = k0;
  m_momenta[m_meshPoints] = m_k0;
}

double NucleonScattering::GetV0(){ return m_V0; }
void NucleonScattering::SetV0(double v){ m_V0 = v; }

double NucleonScattering::GetA(){ return m_a; }
void NucleonScattering::SetA(double a){ m_a = a; }



       /*
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
  // Call Morten's function to fill the weights and k points arrays
  GaussLegendreQuadrature(-1.,1.,m_weights,m_momenta,m_meshPoints);
  // Map w and k from [-1,1] to [0,inf)
  for(int i=0;i<m_meshPoints;i++){
//    printf("Was,\nk = %g\nw = %g\n\n",m_momenta[i],m_weights[i]);
    m_momenta[i] = tan(M_PI_4*(1.+m_momenta[i]));
    double costerm = cos(M_PI_4*(1.+m_momenta[i]));
    m_weights[i] = M_PI_4 * m_weights[i] / (costerm*costerm);
//    printf("Now,\nk = %g\nw = %g\n\n",m_momenta[i],m_weights[i]);
  }

  m_momenta[m_meshPoints] = m_k0;

  m_momentaSet=true;
}


void NucleonScattering::ListPoints(){
  for(int i=0;i<m_meshPoints+1;i++){
    printf("k[%i] = %g\n",i,m_momenta[i]);
  }
}


/* This function sets up the potential matrix for a
 * finite sphere of radius a fm and strength V
*/
void NucleonScattering::FiniteSphere(double V,double a){
m_V0 = V;
m_a  = a;
if(m_momentaSet){
  // allocate memory for potential matrix
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
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

//      printf("V[%i][%i] = %g\n",i,j,m_potential[i][j]);
    }
  }

//  free_matrix((void **) m_potential);
  m_potentialSet=true;
}
else if(m_potentialSet){cout << "Potential already set." << endl;return;}
else{cout << "Run MapPoints() first." << endl;return;}
}


void NucleonScattering::PotentialNP(){
  const double u = 0.7;// [1/fm]
  double Va = -10.463/u;// [MeV]
  double  a =  1.*u;
  a*=a;
  double Vb = -1650.6/u;// [MeV]
  double  b = 4.*u;
  b*=b;
  double Vc =  6484.3/u;// [MeV]
  double  c = 7.*u;
  c*=c;

if(m_momentaSet){
  // allocate memory for potential matrix
  m_potential = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      double potA = log(pow(m_momenta[i]+m_momenta[j],2.) + a);
      potA       -= log(pow(m_momenta[i]-m_momenta[j],2.) + a);
      potA       *= Va / (4. * m_momenta[i]*m_momenta[j]);

      double potB = log(pow(m_momenta[i]+m_momenta[j],2.) + b);
      potB       -= log(pow(m_momenta[i]-m_momenta[j],2.) + b);
      potB       *= Vb / (4. * m_momenta[i]*m_momenta[j]);

      double potC = log(pow(m_momenta[i]+m_momenta[j],2.) + c);
      potC       -= log(pow(m_momenta[i]-m_momenta[j],2.) + c);
      potC       *= Vc / (4. * m_momenta[i]*m_momenta[j]);

      m_potential[i][j] = potA + potB + potC;
      printf("V[%i][%i] = %g\n",i,j,m_potential[i][j]);
    }
  }
  m_potentialSet=true;
}
else if(m_potentialSet){cout << "Potential already set." << endl;return;}
else{cout << "Run MapPoints() first." << endl;return;}
}


void NucleonScattering::SetMatrixA(){
if(m_potentialSet){
  const double TwoOverPi = 2. / M_PI;
  double u0 = 0.;
  // allocate memory for [A]
  m_A = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  for(int i=0;i<m_meshPoints;i++){
    for(int j=0;j<m_meshPoints;j++){
      // Calculate uj terms
      double num   = m_mass* m_weights[j] * m_momenta[j]*m_momenta[j];
      double denom = (m_k0*m_k0 - pow(m_momenta[j],2.));
      double uj = TwoOverPi * num / denom;
      // Calculate elements of [A]
      m_A[i][j] = (-1.)*m_potential[i][j]*uj;
    }// Close j
    m_A[i][i] += 1.;

    // Run sum over momenta for N+1 u term
    u0 += m_k0*m_k0*m_weights[i]/(m_k0*m_k0 - pow(m_momenta[i],2.));
  }// Close i

  u0 *= (-1.)*TwoOverPi*m_mass;
  m_A[m_meshPoints][m_meshPoints] = 1. - m_potential[m_meshPoints][m_meshPoints]*u0;
  m_aSet=true;
}
}


void NucleonScattering::PrintA(){
if(m_aSet){
  for(int i=0;i<=m_meshPoints;i++){
    for(int j=0;j<=m_meshPoints;j++){
      printf("A[%i][%i] = %g\n",i,j,m_A[i][j]);
    }
  }
}
else {cout << "[A] not set." << endl;return;}
}


void NucleonScattering::PhaseShift(){
if(m_aSet){
  // invert [A]
  inverse(m_A,m_meshPoints+1);
  // allocate memory for [R] and initialize
  m_R = (double **) matrix(m_meshPoints+1,m_meshPoints+1,sizeof(double));
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      m_R[i][j] = 0.;
    }
  }
  // perform matrix multiplication
  for(int i=0;i<m_meshPoints+1;i++){
    for(int j=0;j<m_meshPoints+1;j++){
      for(int k=0;k<m_meshPoints+1;k++){
        m_R[i][j] += m_A[i][k] * m_potential[k][j];
      }
    }
  }

  printf("R(k0,k0) = %g\n",m_R[m_meshPoints][m_meshPoints]);
  double delta = (-1.)*m_R[m_meshPoints][m_meshPoints]*m_mass*m_k0;
  delta = atan(delta);
  printf("delta = %g\n",delta);
}
else{cout << "Run SetMatrixA() first." << endl;return;}
}


void NucleonScattering::AnalyticPhaseShift(){
/*
  double K0  = sqrt(2.*m_mass*m_V0)/197.;

  double eta = 0.;
  double delta=0.;

  if(m_k0 < K0){
    cout << "k0 < K" << endl;
    eta = sqrt(K0*K0 - m_k0*m_k0);
  }
  else if(m_k0 > K0){
    cout << "k0 > K" << endl;
    eta = sqrt(m_k0*m_k0 - K0*K0);
  }
  else{cout << "E = V." << endl;}

  double term1 = tanh(eta*m_a);
  double term2 = m_k0 / eta;
  delta = atan(term1*term2) - m_k0*m_a;

  printf("K = %g, a = %g, k = %g\ndelta =% g\n",K0,m_a,m_k0,delta);
*/
  double E = (m_k0*m_k0) / (2. * m_mass);
  double delta = sqrt(E/(E+m_V0));
  delta *= tan(m_a*sqrt(2.*m_mass*(E+m_V0)));
  delta = atan(delta);
  delta -= m_a*sqrt(2.*m_mass*E);

  printf("m = %g\na = %g\nV0 = %g\n",m_mass,m_a,m_V0);
  printf("k = %g\ndelta = %g\n",m_k0,atan(tan(delta)));
  printf("R(k0,k0)= %g\n",-1.*tan(delta)/(m_mass*m_k0));
}

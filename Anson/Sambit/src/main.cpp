#include <dcosmology.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <dnumrecipes.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>


//Cosmology         
const double om0 = 0.27;
const double lam0=0.73;
const double omb = 0.044;
const double hh = 0.701;
const double s8 = 0.8;
const double ns = 0.96;
const double omNu = 0.0;


//Spherical collapse critical (linear) overdensity
const double dcrit=1.68;


//Minimum halo mass for star formation
const double Mmin = 0.8e9; //In solar masses



//Evaluated at delta_l = 0, where delta_l is the linear density field at z=0
//(Note this is different than the convention where delta_l is linearly
//extrapolated to the present day)
//Mmin is the minimum halo mass for star formation
//Here we also assume that the large-scale volume is so large that sigma^2(long)
//is much smaller than sigma^2(Mmin)
double dxHI_ddelta(double Mmin, double z, Cosmology &cc)
{
  double ans = 0.;
  double dc = dcrit;
  double Smin = pow(cc.growthFac(z)*cc.sigma0fM(Mmin,ans,0),2.0);

  ans = -1.*exp(-1.0*dc*dc/2.0/Smin)*sqrt(2.0/M_PI/Smin);
 
  return ans;
}


double fcollapse(double Mmin, double z, Cosmology &cc)
{
  double ans = 0.;
  double dc = dcrit;
  double Smin = pow(cc.growthFac(z)*cc.sigma0fM(Mmin,ans,0),2.0);

  ans = erfc(sqrt(1.0*dc*dc/2.0/Smin));
 
  return ans;
}



using namespace std;

int main(int argc, char *argv[]) {

  Cosmology cc(om0,lam0,omb,hh,s8,ns,omNu);
  double redshift = 9.164, xfrac = 0.02074;

  if (argc>=3) {
     redshift = atof(argv[1]);
     xfrac    = atof(argv[2]);
  }
  cout << dxHI_ddelta(Mmin,redshift,cc) << endl;
  cout << xfrac/fcollapse(Mmin,redshift,cc) << endl;
  cout << fcollapse(Mmin,redshift,cc) << endl;
  cout << xfrac/0.0052 << endl;

  ofstream myfile;
  myfile.open ("info.txt");
  myfile << dxHI_ddelta(Mmin,redshift,cc) << endl;
  myfile << xfrac/fcollapse(Mmin,redshift,cc) << endl;
  myfile.close();

  return 0;

}

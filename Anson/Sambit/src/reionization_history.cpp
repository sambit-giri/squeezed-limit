#include <dcosmology.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <dnumrecipes.h>



//Cosmology         
const double om0 = 0.28;
const double lam0=0.72;
const double omb = 0.046;
const double hh = 0.7;
const double s8 = 0.817;
const double ns = 0.96;
const double omNu = 0.0;


//Spherical collapse critical (linear) overdensity
const double dcrit=1.68;


//Minimum halo mass for star formation
const double Mmin = 2.e9; //In solar masses



//Evaluated at delta_l = 0, where delta_l is the linear density field at z=0
//(Note this is different than the convention where delta_l is linearly
//extrapolated to the present day)
//Mmin is the minimum halo mass for star formation
//Here we also assume that the large-scale volume is so large that sigma^2(long)
//is much smaller than sigma^2(Mmin)

double xHI_delta(double Mmin, double z, Cosmology &cc)
{
  double ans = 0.;
  double dc = dcrit;
  double Smin = pow(cc.growthFac(z)*cc.sigma0fM(Mmin,ans,0),2.0);
  double efficiency = 1.0;

  ans = sqrt(1.0*dc*dc/2.0/Smin);
  ans = erf(ans);
  ans *= efficiency;
 
  return ans;
}




using namespace std;

int main(int argc, char *argv[]) {


  Cosmology cc(om0,lam0,omb,hh,s8,ns,omNu);
  double redshift_i = 6.;
  int i;

  for(i=0;i<15;i++) {
    cout << xHI_delta(Mmin,redshift_i+i,cc) << endl;
  }

  return 0;

}

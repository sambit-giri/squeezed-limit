// dcosmology.h
// Includes constants describing the overall cosmology of a universe
// and prototypes for class Cosmology

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

// The following utility constants are in cgs units.
const double CRITDENSITY(1.8791e-29);   // Current critical density, 
                                        // in g/cm^3, not including h^2!
const double CRITDENMSOLMPC(2.7755e11); // Current critical density, 
                                        // in Msol/Mpc^3, not including h^2!
const double MPROTON(1.6726e-24);           // in grams
const double BOLTZK(1.3806e-16);            // in erg/K
const double EDDINGTONLUMINOSITY(87.83471);   // For 1 Msol, in erg/s; enter 
                                             // ln of 1.4e38
const double YEAR(3.156e7);                 // # seconds/year
const double KPC(3.09e21);                  // cm/kpc
const double UNH(3.241e-18);                // In 1/s
const double TCMB(2.726);                   // CMB temperature at z=0
const double MUB(1.22);           // Mean molecular weight, for primordial
const double ELECTRONCHARGE(4.80e-10);
const double MELECTRON(9.11e-28);
const double LIGHTSPEED(3.0e10);

// Constant needed for Press-Schechter calculations
const int MAXBINS(7);

class Cosmology {

 public:
  Cosmology(double om0, double lam0, double omb, double h, double s8,
	    double n, double omNu);
  ~Cosmology();

  // General cosmology functions
  double omegaZ(double zCurrent);
  double lambdaZ(double zCurrent);
  double hubbleZ(double zCurrent);
  double rhoCritZ(double zCurrent);
  double drdz(double zCurrent);
  double cosmicTime(double zCurrent);
  double zFromT(double ht0, double zGuess=100);
  double dzdh0t(double h0time);
  double coordDistance(double zCall);
  double confTime(double zCall);
  double lumDistance(double zCall);
  double angDiamDistance(double zCall);
  double nh(double zCall);

  // Perturbation theory functions
  double DeltaC(double zCurrent);
  double delCrit0(double z);
  double growthFac(double z);

  // Finding mass scale corresponding to sigma
  double findMassAtSigma(double n, double z0, double deltax = -1.0);
  
  // Finding collapse fraction
  double fCollPSExact(double z, double mMin = -1.0);
  double fColl(double z, double mMin = -1.0, int massFcn = 0);
  double nCollObject(double z, double mMin = -1.0);

  // Linear bias
  double biasPS(double z, double mass);

  // Press-Schechter functions and extensions
  double dndlM(double z, double tM);
  double dndlMSheth(double z, double tM);
  double dndlMJenkins(double z, double tM);
  void resetPowerSpectrum();
  double sigma0fM(double tM, double &dsdM, int iDeriv);

  // Power spectrum evaluation
  double powerSpectrum(double k);
  void TFSetParameters(void);
  double TFMaster(double k, double &dTFdk,int iDeriv);
  double TF_BBKS(double k);

  // Extended Press-Schechter functions
  double mergerRate(double mTot, double mAccrete, double z);
  double mergerKernelSym(double m1, double m2, double z);
  double mergerKernel(double m1, double m2, double z);
  double deltacderiv(double z);
  double formTimeCDF(double zform, double mass, double zfinal);
  double formTimePDF(double zform, double mass, double zfinal);
  double medianFormTime(double mass, double zfinal);
  double collapseRedshift(double mass, double zfin, double fInput);

  // Data member retrieval functions
  double getOmega0();
  double getOmegam();
  double getOmegab();
  double getOmbhh();
  double getOm0hh();
  double getLambda0();
  double getH();
  double getnSpec();
  double getScale();
  double getShape();
  double getSnorm();

 protected:

  // Parameters of universe
  double omega0;
  double lambda0;                   // Actually omega(sub)lambda
  double h;
  double omegab;
  double ombhh;
  double omegam;
  double om0hh;
  double nspec;                    // Power spectrum index
  double shape;                    // Gamma
  double omeganu;
  double sig8;
  // Some more variables, used and calculated by Press-Schechter functions
  // but globally useful as constants.
  double thetaCMB;
  double scale;
  double sNorm;
  double soundHorizon;
  double zEquality;
  double alphaNu;
  double betaC;
};

// Finding mass scale corresponding to sigma
double setMassSigma(double mass, double n, double z, double deltax, 
		   Cosmology *c1, int flag);
double massSigma(double mass);

// Finding collapse fraction
double setFCollInt(double mass, Cosmology *c1, double z1, int massFcn, 
		  int flag);
double fCollInt(double mass);
double setNCollInt(double mass, Cosmology *c1, double z1, int flag);
double nCollInt(double mass);

// Integrands for calculating PS functions
double sigmatop2(double k);
double dsigmatop2(double k);
double sigmatop(double kl);
double setSigmatop(double kl, Cosmology *c1, int flag);
double dsigmatop(double kl);
double setdsigmatop(double kl, Cosmology *c1, int flag);

// Derivative functions for EPS
double setDeltacritDeriv(double z, Cosmology *cos, int flag);
double deltaCritDeriv(double z);

// Functions for EPS formation time calculation
void setFormTimeCDFInt(double stilde, double y[], double result[], 
		       Cosmology *cos, double omt, double m, 
		       double shalf, double stwo, int flag);
void formTimeCDFInt(double s, double y[], double result[]);
double setInvertSigmaM(double mass, Cosmology *cos, double sig, int flag);
double invertSigmaM(double mass);
void setFormTimePDFInt(double stilde, double y[], double result[], 
		       Cosmology *cos, double omt, double m, 
		       double shalf, double stwo, int flag);
void formTimePDFInt(double s, double y[], double result[]);
double setFindMedian(double zform, double m, double zfin,
		     Cosmology *cos, int flag);
double findMedian(double zform);

/* Below are some ancillary functions that are also generally useful */

/* Cosmology-dependent halo characteristics. */
double jeansMass(Cosmology *c, double z);
double filterMass(Cosmology *c, double z);
double minIonMass(Cosmology *c, double z);
double coolMass(Cosmology *c, double z);

double rvir(Cosmology *c, double mass, double z);
double vcirc(Cosmology *c, double mass, double z);
double tvir(Cosmology *c, double mass, double z, double mu);
double MfromTvir(Cosmology *c, double tv, double z, double mu);
double RComfromM(double m, Cosmology *c);
double MfromRCom(double R, Cosmology *c);

/* Interpolation functions.  Need to set up sigm if going to be used! */
double nm(double tM, double z, Cosmology *c);
double fcoll(double z, double mMin, Cosmology *c);
double nmcond(double tM, double z, double mBubble, double deltaBubble, 
	      Cosmology *c);
double nm_st(double tM, double z, Cosmology *c);
double nm_st_fit(double tM, double z, Cosmology *c);
double nmcond_st(double tM, double z, double mBubble, double deltaBubble, 
		 Cosmology *c);
double biasm(double mass, double z, Cosmology *c);
double biasm_st(double mass, double z, Cosmology *c);
double sigm(double m, double &dsdm, Cosmology *c, int flag);

#endif


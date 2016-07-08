#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;

/*======================================================================================================================== Emperical Neutron Star parameters - working in units of MeV for energies and fm for distances ========================================================================================================================*/

const double MeVtoinvFM = 0.00507;				// approximate conversion factor
const double rho0 = 0.153;					// staturation density (fm^-3)
const double kFermi = 1.31;					// Fermi momentum (fm^-1)
const double kFermi2 = 1.31*1.31;				// Fermi momentum (squared)
const double Kompress = 200*MeVtoinvFM;				// Compression modulus (fm^-1)
const double mass = 4.758;					// nucleon mass (fm^-1) 
const double mstar = 0.7*mass;					// effective mass
const double asymm = 32.5*MeVtoinvFM;		// symmetry energy coefficient (fm^-1)
const double BperA = -16.3*MeVtoinvFM;				// Binding energy (fm^-1)
const double pi = M_PI;						// Pie (yummm)

struct arg_params{
	double mass; 
};

/*================================================================================================================================================= Constant parameters at saturation ===========================================================================================================================================================*/
double Gomega2(double RHO){
	double fieldstrength = (mass + BperA - sqrt(kFermi2 + mstar * mstar))/RHO;	
	return fieldstrength;
}
double alpha1(double RHO){
	double a1 = Kompress - Gomega2(RHO)*6*pow(kFermi, 3)/(pi*pi) - 3 * kFermi2/sqrt(kFermi2 + mstar*mstar);
	return a1;
}
double beta1(double RHO){
	double b1 = 2*pow((mass - mstar),2)*alpha1(RHO);	//NOTE: I CHANGED WHAT beta1 IS FROM LIT.
	return b1;
}
double gamma1(double RHO){
	double g1 = 3*pow((mass - mstar),2) * alpha1(RHO);
	return g1;
}

double alpha2 = 0.5 * pow((mass - mstar),2);
double beta2 = (1/3) * mass * pow((mass - mstar),3);
double gamma2 = 0.25 * pow((mass - mstar),4);

double alpha3 = (mass - mstar);
double beta3 = mass * pow((mass - mstar),2);
double gamma3 = pow((mass - mstar),3);

double Xfac = kFermi/mstar;
double Tfac = sqrt(1 + Xfac*Xfac);

double I1 = 2/(pi*pi) * mstar*mstar* (0.5*Xfac*Tfac + Xfac/Tfac - 1.5*log(Xfac+Tfac));
double I2 = 1/(2*pi*pi) * pow(mstar, 4) * (Xfac * pow(Tfac, 3) - 0.5 * Xfac * Tfac - 0.5 * log(Xfac + Tfac));
double I3 = 1/(pi*pi) * pow(mstar, 3) * ( Xfac * Tfac - log(Xfac + Tfac));

double del1(double RHO){
	double d1 = -alpha1(RHO)*I1 - 6 * kFermi2*kFermi/(pi*pi) * (mstar*mstar/(kFermi2 + mstar*mstar));
	return d1;
}

double del3 = I3;
double del2(double RHO){
	double d2 = RHO* (mass + BperA) - I2 - 0.5 * Gomega2(RHO) * pow(RHO,2);
	return d2;
}


/*===================================================================================================================================================== variable parameters =====================================================================================================================================================================*/

double C(double RHO){
	double b = (alpha2*del1(RHO) - alpha1(RHO)*del2(RHO));
	b *= 1/(alpha2*gamma1(RHO) - alpha1(RHO)*gamma2);
	b *= (alpha3*del1(RHO) - alpha1(RHO)*del3);
	b *=1/(alpha3*gamma1(RHO) - alpha1(RHO)*gamma3);
	return b;
}

double B(double RHO){
	double c = (alpha2*del1(RHO) - alpha1(RHO)*del2(RHO))-(alpha2*gamma1(RHO) - alpha1(RHO)*gamma2);
	c *= 1/(alpha2*beta1(RHO) - alpha1(RHO)*beta2);
	c *= C(RHO);
	return c;
}

double Grho2(double RHO){
	double fieldstrength = asymm - kFermi2/(6*sqrt(kFermi2 + mstar*mstar));
	fieldstrength *= 8/RHO;
	return fieldstrength;
}

double Gsigma2(double RHO){
	double fieldstrength = alpha1(RHO)/(del1(RHO) - gamma1(RHO)*C(RHO) - beta1(RHO)*B(RHO));
	return fieldstrength;
}

/*================================================================================================================================================ functions to be integrated ===================================================================================================================================================================*/

double rhoS(double kF, void *param){				//Scalar density Integrand
	double Integrand = 0;
	struct arg_params *fp = (struct arg_params*) param;
	double mst = fp->mass;

	double kF2 = kF*kF;
	double kfactor = sqrt(kF2 + mst*mst);

	Integrand = 2*kF2*mst/(pi*pi * kfactor);
	return Integrand;
}

double goko(double kF, void *param){				//Kinetic Energy Integrand
	double Integrand = 0;
	struct arg_params *fp = (struct arg_params*) param;
	double mst = fp->mass;

	double kF2 = kF*kF;
	double kfactor = sqrt(kF2 + mst*mst);

	Integrand = (2/(pi*pi)) * kF2* kfactor;
	return Integrand;
}


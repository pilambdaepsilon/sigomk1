#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;

/*===============================================================================================================
========= Emperical Nuclear parameters - working in units of MeV for energies and fm for distances ==============
===============================================================================================================*/
const double MeVtoinvFM = 0.00507;				// approximate conversion factor
const double rho0 = 0.153;					// staturation density (fm^-3)
const double kFermi = 1.31;					// Fermi momentum (fm^-1)
const double kFermi2 = pow(kFermi, 2);				// Fermi momentum (squared)
const double Kompress = 200*MeVtoinvFM;				// Compression modulus (fm^-1)
double mass = 938*MeVtoinvFM;					// nucleon mass (fm^-1) 
double mstar = 0.7*mass;					// effective mass
const double asymm = 32.5*MeVtoinvFM;				// symmetry energy coefficient (fm^-1)
const double BperA = -16.3*MeVtoinvFM;				// Binding energy (fm^-1)
const double pi = M_PI;						// Pie (yummm)


/*===============================================================================================================
================================== Constant parameters at saturation ============================================
===============================================================================================================*/
double Gomega2(double RHO){
	double fieldstrength = (mass + BperA - sqrt(kFermi2 + mstar * mstar))/RHO;	
	return fieldstrength;
}

double alpha1(double RHO){
	double a1 = Kompress - Gomega2(RHO)*6*pow(kFermi, 3)/(pi*pi) - 3 * kFermi2/sqrt(kFermi2 + mstar*mstar);
	return a1;
}
double beta1(double RHO){
	double b1 = 2*mass*(mass - mstar)*alpha1(RHO);	//NOTE: I CHANGED WHAT beta1 IS FROM LIT.
	return b1;
}
double gamma1(double RHO){
	double g1 = 3*pow((mass - mstar),2) * alpha1(RHO);
	return g1;
}
double del1(double RHO, double IONE){
	double d1 = -alpha1(RHO)*IONE - 6 * kFermi2*kFermi/(pi*pi) * (mstar*mstar/(kFermi2 + mstar*mstar));
	return d1;
}


double alpha2 = (mass - mstar);
double beta2 = mass * pow((mass - mstar),2);
double gamma2 = pow((mass - mstar),3);
double del2(double ITWO){
	double d2 = ITWO;
	return ITWO;
}

double alpha3 = 0.5 * pow((mass - mstar),2);
double beta3 = (1/3) * mass * pow((mass - mstar),3);
double gamma3 = 0.25 * pow((mass - mstar),4);
double del3(double RHO, double ITHREE){
	double d3 = RHO* (mass + BperA) - ITHREE - 0.5 * Gomega2(RHO) * pow(RHO,2);
	return d3;
}



double Xfac = kFermi/mstar;
double Tfac = sqrt(1 + Xfac*Xfac);
double ithree = 1/(2*pi*pi) * pow(mstar, 4) * (Xfac * pow(Tfac, 3) - 0.5 * Xfac * Tfac - 0.5 * log(Xfac + Tfac));

/*===============================================================================================================
====================================== variable parameters ======================================================
===============================================================================================================*/
double B2(double RHO, double IONE, double ITWO, double ITHREE){
	double ag21 = ( alpha2*gamma1(RHO) - alpha1(RHO)*gamma2 );
	double ag32 = ( alpha3*gamma2 - alpha2*gamma3 );
	
	double ad21 = ( alpha2*del1(RHO, IONE) - alpha1(RHO)*del2(ITWO) );
	double ad32 = ( alpha3*del2(ITWO) - alpha2*del3(RHO, ITHREE) );

	double ab21 = ( alpha2*beta1(RHO) - alpha1(RHO)*beta2 );
	double ab32 = ( alpha3*beta2 - alpha2*beta3 );

	double bnumerator = ag21*ad32 - ag32*ad21;
	double bdenominator = ag21*ab32 - ag32*ab21;

	double b = bnumerator/bdenominator;
	return b;
}

double C2(double RHO, double IONE, double ITWO, double ITHREE){
	double ag21 = ( alpha2*gamma1(RHO) - alpha1(RHO)*gamma2 );
	double ag32 = ( alpha3*gamma2 - alpha2*gamma3 );
	
	double ad21 = ( alpha2*del1(RHO, IONE) - alpha1(RHO)*del2(ITWO) );
	double ad32 = ( alpha3*del2(ITWO) - alpha2*del3(RHO, ITHREE) );

	double ab21 = ( alpha2*beta1(RHO) - alpha1(RHO)*beta2 );
	double ab32 = ( alpha3*beta2 - alpha2*beta3 );

	double cnumerator = ad32 - ab32*B2(RHO, IONE, ITWO, ITHREE);
	double cdenominator = ag32;

	double c = cnumerator/cdenominator;
	return c;
}


double C(double RHO, double IONE, double ITWO, double ITHREE){
	double ag21 = ( alpha2*gamma1(RHO) - alpha1(RHO)*gamma2 );
	double ag31 = ( alpha3*gamma1(RHO) - alpha1(RHO)*gamma3 );
	
	double ad31 = ( alpha3*del1(RHO, IONE) - alpha1(RHO)*del3(RHO, ITHREE) );
	double ad21 = ( alpha2*del1(RHO, IONE) - alpha1(RHO)*del2(ITWO) );

	double ab31 = ( alpha3*beta1(RHO) - alpha1(RHO)*beta3 );
	double ab21 = ( alpha2*beta1(RHO) - alpha1(RHO)*beta2 );

	double cnumerator = ab31*ad21*ab21*ad31;
	double cdenominator = ab31*ag21*ab21*ag31;

	double c = cnumerator/cdenominator;
	return c;
}

double B(double RHO, double IONE, double ITWO, double ITHREE){
	double ag21 = ( alpha2*gamma1(RHO) - alpha1(RHO)*gamma2 );
	double ag31 = ( alpha3*gamma1(RHO) - alpha1(RHO)*gamma3 );
	
	double ad31 = ( alpha3*del1(RHO, IONE) - alpha1(RHO)*del3(RHO, ITHREE) );
	double ad21 = ( alpha2*del1(RHO, IONE) - alpha1(RHO)*del2(ITWO) );

	double ab31 = ( alpha3*beta1(RHO) - alpha1(RHO)*beta3 );
	double ab21 = ( alpha2*beta1(RHO) - alpha1(RHO)*beta2 );

	double bnumerator = (ad21 - ag21)*C(RHO, IONE, ITWO, ITHREE);
	double bdenominator = ab21;

	double b = bnumerator/bdenominator;
	return b;
}

/*===================================================================================================*/
double Grho2(double RHO){
	double fieldstrength = asymm - kFermi2/(6*sqrt(kFermi2 + mstar*mstar));
	fieldstrength *= 8/RHO;
	return fieldstrength;
}

double Gsigma2(double RHO, double CONSTB, double CONSTC, double IONE){
	double sigmanumerator = alpha1(RHO);
	double sigmadenominator = del1(RHO, IONE) - beta1(RHO)*CONSTB - gamma1(RHO)*CONSTC;

	double fieldstrength = sigmanumerator/sigmadenominator;
	return fieldstrength;
}

/*===============================================================================================================
================================= functions to be integrated ====================================================
===============================================================================================================*/
double rhoS(double k, void *param){			//Scalar density Integrand
	double Integrand = 0;
	double m = *(double *) param;

	double k2 = k*k;
	double E = sqrt(k2 + m*m);

	Integrand = 2*k2*m/(pi*pi * E);
	return Integrand;
}

double I1(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double E = sqrt(k*k + m*m);			

	Integrand = (2/(pi*pi)) * pow(k,4) / pow(E,3);
	return Integrand;
}

double I2(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = (2/(pi*pi)) * k2* m / E;
	return Integrand;
}

double goko(double k, void *param){			//Kinetic Energy Integrand-same form as I3
							//so use it in the integratori to solve for del3
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);		

	Integrand = (2/(pi*pi)) * k2* E;
	return Integrand;
}


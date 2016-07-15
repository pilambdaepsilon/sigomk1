#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
using namespace std;

/*==============================================================================================================
==== Emperical Nuclear parameters - working in units of MeV and fm^-1 for energies and fm for distances ========
==============================================================================================================*/
const double pi = M_PI;						// Pie (yummm)
const double MeVtoinvFM = 0.0008065*(2*pi);			// approximate conversion factor
const double rho0 = 0.153;					// staturation density (fm^-3)
double Kompress = 250*MeVtoinvFM;				// Compression modulus (fm^-1)
double mass = 938*MeVtoinvFM;					// nucleon mass (fm^-1) 
double mstar = 0.75*mass;					// effective mass
const double asymm = 32.5*MeVtoinvFM;				// symmetry energy coefficient (fm^-1)
const double BperA = -16.3*MeVtoinvFM;				// Binding energy (fm^-1)

/*===============================================================================================================
================================== Constant parameters at saturation ============================================
===============================================================================================================*/
double Gomega2(double RHO, double kF){
	double fieldstrength = (mass + BperA - sqrt(kF*kF + mstar * mstar))/RHO;	
	return fieldstrength;
}
double Grho2(double RHO, double kF){
	double fieldstrength = 8*(asymm - kF*kF/(6*sqrt(kF*kF + mstar*mstar)))/RHO;
	return fieldstrength;
}

double alpha1(double RHO, double kF){
	double a1 = Kompress - Gomega2(RHO, kF)*6*pow(kF, 3)/(pi*pi) - 3 * kF*kF/sqrt(kF*kF + mstar*mstar);
	return a1;
}
double beta1(double RHO, double kF){
	double b1 = 2*mass*pow((mass - mstar),2)*alpha1(RHO,kF);	
	return b1;
}
double gamma1(double RHO, double kF){
	double g1 = 3*pow((mass - mstar),2) * alpha1(RHO, kF);
	return g1;
}
double del1(double RHO, double IONE, double kF){
	double d1 = -alpha1(RHO, kF)*IONE - 6 * pow(kF,3)/(pi*pi) * (mstar*mstar/(kF*kF + mstar*mstar));
	return d1;
}


double alpha3 = (mass - mstar);
double beta3 = mass * pow((mass - mstar),2);
double gamma3 = pow((mass - mstar),3);
double del3(double ITHREE){
	double d3 = ITHREE;
	return ITHREE;
}

double alpha2 = 0.5 * pow((mass - mstar),2);
double beta2 = (1/3) * mass * pow((mass - mstar),3);
double gamma2 = 0.25 * pow((mass - mstar),4);
double del2(double RHO, double ITWO, double kF){
	double d2 = RHO* (mass + BperA) - ITWO - 0.5 * Gomega2(RHO, kF) * pow(RHO,2);
	return d2;
}



/*===============================================================================================================
====================================== variable parameters ======================================================
===============================================================================================================*/

double C(double RHO, double IONE, double ITWO, double ITHREE, double kF){
	double ag21 = ( alpha2*gamma1(RHO, kF) - alpha1(RHO, kF)*gamma2 );
	double ag31 = ( alpha3*gamma1(RHO, kF) - alpha1(RHO, kF)*gamma3 );
	
	double ad31 = ( alpha3*del1(RHO, IONE, kF) - alpha1(RHO, kF)*del3(ITHREE) );
	double ad21 = ( alpha2*del1(RHO, IONE, kF) - alpha1(RHO, kF)*del2(RHO, ITWO, kF) );

	double ab31 = ( alpha3*beta1(RHO, kF) - alpha1(RHO, kF)*beta3 );
	double ab21 = ( alpha2*beta1(RHO, kF) - alpha1(RHO, kF)*beta2 );

	double cnumerator = ab31*ad21*ab21*ad31;
	double cdenominator = ab31*ag21*ab21*ag31;

	double c = cnumerator/cdenominator;
	return c;
}

double B(double RHO, double IONE, double ITWO, double ITHREE, double kF){
	double ag21 = ( alpha2*gamma1(RHO, kF) - alpha1(RHO, kF)*gamma2 );
	double ag31 = ( alpha3*gamma1(RHO, kF) - alpha1(RHO, kF)*gamma3 );
	
	double ad31 = ( alpha3*del1(RHO, IONE, kF) - alpha1(RHO, kF)*del3(ITHREE) );
	double ad21 = ( alpha2*del1(RHO, IONE, kF) - alpha1(RHO, kF)*del2(RHO, ITWO, kF) );

	double ab31 = ( alpha3*beta1(RHO, kF) - alpha1(RHO, kF)*beta3 );
	double ab21 = ( alpha2*beta1(RHO, kF) - alpha1(RHO, kF)*beta2 );

	double bnumerator = (ad21 - ag21)*C(RHO, IONE, ITWO, ITHREE, kF);
	double bdenominator = ab21;

	double b = bnumerator/bdenominator;
	return b;
}


double Gsigma2(double RHO, double CONSTB, double CONSTC, double IONE, double kF){
	double sigmanumerator = alpha1(RHO, kF);
	double sigmadenominator = del1(RHO, IONE, kF) - beta1(RHO, kF)*CONSTB - gamma1(RHO, kF)*CONSTC;

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
double eye1(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
	return 2/(pi*pi) * m*m*(0.5*x*t + x/t - (1.5)*log(x+t));
}
	 
double eye2(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
	return 1./(pi*pi) * pow(m, 3)*(x*t - log(x+t));
}
double eye3(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
	return 1./(2*pi*pi) * pow(m, 4)*(-0.5*x*t + x*pow(t, 3) - 0.5*log(x+t));
}
	


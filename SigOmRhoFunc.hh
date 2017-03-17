#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"
#include <vector>
using namespace std;

/*==============================================================================================================
==== Emperical Nuclear parameters - working in units of MeV and fm^-1 for energies and fm for distances ========
==============================================================================================================*/
const double MeVtoinvFM = 0.0008065*(2*pi);			// approximate conversion factor
//const double rho0 = 0.153;					// staturation density (fm^-3)
//double Kompress = 200*MeVtoinvFM;				// Compression modulus (fm^-1)
//double mass = 938*MeVtoinvFM;					// nucleon mass (fm^-1) 
//double massOmega = 783*MeVtoinvFM;				// omega mass (fm^-1) 
//double massSigma = 550*MeVtoinvFM;				// approximate sigma mass (fm^-1) 
//double mstar = 0.75*mass;					// effective mass
//const double asymm = 32.5*MeVtoinvFM;				// symmetry energy coefficient (fm^-1)
//const double BperA = -16.3*MeVtoinvFM;				// Binding energy (fm^-1)
const double MPI = 90e3*MeVtoinvFM;
/*===============================================================================================================
================================= functions to be integrated ====================================================
===============================================================================================================*/
double I3(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double E = sqrt(k*k + m*m);			

	Integrand = pow(k,4) / pow(E,3);
	return Integrand;
}

double I2(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* E;
	return Integrand;
}

double I1(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* m / E;
	return Integrand;
}

double IP(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* k2 / E;
	return Integrand;
}

/*=================================== DM CONSTRAINED PARAMETERS ===========================================*/
double Gchi(double sigchi2, double mchi, double mpi){
	double arg = sqrt(sigchi2) * mpi*mpi/( sqrt(20.)* mchi);
	return sqrt(arg);
}
double Gpi(double sigchin, double mn, double mchi, double mpi, double gchi){
	double arg = sqrt(sigchin)* (mn + mchi)*mpi*mpi/( sqrt(3.)*gchi*mn*mchi);
	return sqrt(arg);
}


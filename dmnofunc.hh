//#include <gsl/gsl_integration.h>
//#include "LHAPDF/LHAPDF.h"
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include "CONSTANTSandCONVERSIONS.hh"

/*====================== EMPERICAL NS PARAMETERS =======================*/
const double Msolar = 1.98892e30*convkgtoGeV;						// mass of Sol in GeV
double TNS = 1e5*convKtoGeV;								// NS temperature in K->GeV
double TDM = TNS;
double RNS = 11.914e5*convcmtoinvGeV;							// NS radius in cm->GeV^-1
double MNS = 2.39*Msolar;								// NS mass in GeV
double vbar = 220e5*convcmtoinvGeV/convstoinvGeV;					// mean thermal DM speed in cm/s->fraction of c - MB distrib
double v2 = vbar*vbar;									// the square pops up more frequently
const double Nmass = 0.94;								// Neutron mass in GeV
const double rhoDM = 1./pow(convcmtoinvGeV,3.);						// GeV/cm^3 -> GeV^4 - mean DM density at 1 kpc from Galactic center
const double rho_baryons = 2.4911e12*convkgtoGeV/pow(convcmtoinvGeV, 3.);			// mean NS baryon mass density in kg/cm^3->GeV^4
const double CoreDens = rho_baryons/(Nmass);						// mean NS baryon number density in GeV^3

//double pFermi = Planckbar*pow((3*pi*pi*CoreDens), 1./3.)/clight;//GeV
double pFermi = 0.426;									// Fermi momentum, neutrons, in GeV
double T1 = 3.5;									// containment time	
double TMAX = 1e10;


double tcon(double mDM, double sigchin){
	double sigcrit = Nmass*RNS*RNS/MNS;
	return 3.*pi*mDM*pow(RNS, 3./2.)*sigcrit/(4.*Nmass*sqrt(2.*GNewton0*MNS)*sigchin*pow(convcmtoinvGeV,2.))*sqrt(mDM/Nmass)/convstoinvGeV/(60*60*24*365);
}
double radius(double mDM, double mred, double t, double sigchin){
	double eta = 8*pi*sqrt(2.) * pow(mred, 3.) * pow(CoreDens, 2.) * sigchin*pow(convcmtoinvGeV,2.) *RNS*RNS* GNewton0 / (3.*mDM * pFermi);
	double r = RNS/sqrt(1.+eta*(t-T1)*365*24*60*60*convstoinvGeV);
	return r;
}	

/*============= ORIGINAL EVOLUTION =============*/
double Nchi0(double Cc, double Cs, double time){
	return Cc/Cs * (exp(Cs * time) - 1.);
}

/* ============ GEOMETRIC ======================*/
double tgeo(double mDM, double mred, double sigchin, double sigchi2, double Cc, double Cs){ 
	double time = T1;
	double tGEO = TMAX*2;
	while(time <= 1e10){
		double rad = radius(mDM, mred, time, sigchin)/RNS;
		double rfid = sqrt( Nchi0(Cc, Cs, time)* sigchi2*pow(convcmtoinvGeV,2.)/(pi * RNS *RNS));
		if (rad <= 1.001*rfid && rad >= 0.999 *rfid ){
			tGEO = time;	
		}
		time *= 1.0001;
	}
	if(tGEO < TMAX){
		return tGEO;
	}
	else if (tGEO > TMAX){
		return TMAX*2;
	}
}
double Nchigeo(double mDM, double mred, double geotime, double time, double sigchin, double sigchi2, double Cc, double Cs){
	double eta = 8*pi*sqrt(2.) * pow(mred, 3.) * pow(CoreDens, 2.) * sigchin*pow(convcmtoinvGeV,2.) *RNS*RNS* GNewton0 / (3*mDM * pFermi);		//NATURAL
	double radgeo = radius(mDM, mred, geotime, sigchin);												//NATURAL
	double Ngeo = pi* pow(radgeo, 2)/(sigchi2*pow(convcmtoinvGeV,2.));
	return Cc*time + Cs/(60*60*24*365)*convinvstoGeV*Ngeo*RNS*RNS/(pow(radgeo, 2)*eta) * log(1 + eta*(time - T1)*365*24*60*60*convstoinvGeV);
}


/* ============ THERMALIZATION ===================*/
double ttherm(double mDM, double mn, double mred, double sigchin){
	return pFermi/(6.*sqrt(2.)*TNS*CoreDens*sigchin*pow(convcmtoinvGeV,2.)) * pow( (mn/mred), 3)*pow(mDM/mn,2.)/convstoinvGeV/(60*60*24*365);
}
double Nchitherm(double mDM, double mred, double thermtime, double time, double sigchin, double sigchi2, double Cc, double Cs){
	double radtherm = sqrt(9.* TNS/(4.*pi*GNewton0*CoreDens*Nmass*mDM));
//	double radtherm = radius(mDM, mred, thermtime, sigchin);
	return (Cc + Cs*pi*pow(radtherm, 2)/(sigchi2*pow(convcmtoinvGeV,2.) ) )* time;
}
double PauliBlocking(double mred){
	double xi = 2.* sqrt(GNewton*MNS/RNS)*mred;
	if (xi < 1.){
		return xi;}
	else{ return 1.;}
}



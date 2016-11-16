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

using namespace std;

/*======================= UNIT CONVERSION ======================*/
double pi = M_PI;
const double convGeVtokg = 1.78e-27;
const double convGeVtoinvm = 5.08e15;
const double convGeVtoinvs = 1.52e24;
const double convmtoinvGeV = 5.08e15;
const double convstoinvGeV = 1.52e24;
const double convGeVtoK =1.16e13;
const double convkgtoGeV = 5.62e26;
const double convinvmtoGeV = 1/5.08e15;
const double convinvstoGeV = 1/1.52e24;
const double convKtoGeV =1/1.16e13;
const double convGeVtoinvFM = 0.0008065*(2*pi)*1e-3;			

/*===================== CONSTANTS ==========================*/
double Gnewton = 6.67408e-5/(convkgtoGeV);			// Gravitational constant in cm^3/kg/s^2 -> cm^3/GeV/s^2
double kBoltzmann = 1.3806e-19*convkgtoGeV/convKtoGeV;		// cm^2*kg/s^2/K -> cm^2/s^2
double clight = 3e10;						// cm/s
const double Planckbar = 6.626e-30*convkgtoGeV;			// m^2*kg/s -> cm^2*GeV/s

/*====================== EMPERICAL NS PARAMETERS =======================*/
const double Msolar = 1.98e30*convkgtoGeV;	// mass of Sol in kg->GeV
const double rho_baryons = 5.7e11*convkgtoGeV;	// mean NS baryon mass density in kg/cm^3->GeV/cm^3
double TNS = 1e5*convKtoGeV;					// NS temperature in K->GeV
double TDM = TNS;
double RNS = 10.6e5;						// NS radius in cm
double MNS = 1.44*Msolar;					// NS mass in GeV
double vbar = 220e5;						// mean thermal DM speed in cm/s - MB distrib
double v2 = vbar*vbar;						// the square pops up more frequently
const double Nmass = 0.94;					//GeV
const double rhoDM = 1;						//GeV/cm^3 - mean DM density at 1 kpc from Galactic center
const double CoreDens = rho_baryons/(Nmass);

//double pFermi = Planckbar*pow((3*pi*pi*CoreDens), 1./3.)/clight;//GeV
double pFermi = sqrt(2*Nmass*.097);
double T1 = 3.5;						//containment time	
double TMAX = 1e10;

double radius(double mDM, double mred, double t, double sigchin){
	double eta = 8*pi*sqrt(2.) * pow(mred, 3.) * pow(CoreDens, 2.) * sigchin *RNS*RNS* Gnewton / (3*mDM * pFermi) * pow(clight, 0) * pow(Planckbar, 0);
	double r = RNS/sqrt(1+eta*(t-T1));
	return r;
}	
double radius2(double mDM, double mn, double t, double sigchin){
	double eta = mDM/mn * 2.8e10/sigchin;
	return sqrt(eta/t);
}	
double radius3(double mDM, double mn, double t, double sigchin){
	double eta = mn*mn/(mDM*mDM) * 2.8e10/sigchin;
	return sqrt(eta/t);
}	

double Nchi0(double Cc, double Cs, double time){
	return Cc/Cs * (exp(Cs * time) - 1.);
}

/* ============ GEOMETRIC ======================*/
double tgeo(double mDM, double mred, double sigchin, double sigchi2, double Cc, double Cs){ 
	double time = T1;
	double tGEO = 0;
	while(time <= 1e10){
		double rad = radius(mDM, mred, time, sigchin)/RNS;
		double rfid = sqrt( Nchi0(Cc, Cs, time)* sigchi2/(pi * RNS *RNS));
		if (rad <= 1.1*rfid && rad >= 0.9 *rfid ){
			tGEO = time;	
		}
		time *= 1.001;
	}
	return tGEO;
}
double Nchigeo(double mDM, double mred, double geotime, double time, double sigchin, double sigchi2, double Cc, double Cs){
	double radgeo = radius(mDM, mred, geotime, sigchin);
	double Ngeo = pi* pow(radgeo, 2)/sigchi2;
	return Cc*time + Cs*Ngeo*geotime * log10(time);
}


/* ============ THERMALIZATION ===================*/
double ttherm(double mDM, double mn, double mred, double sigchin){
	return mDM*mDM*pFermi/(6*sqrt(2)*TNS*mn*mn*CoreDens*sigchin) * pow( (mn/mred), 3);
}
double Nchitherm(double mDM, double mred, double thermtime, double time, double sigchin, double sigchi2, double Cc, double Cs){
	double radtherm = radius(mDM, mred, thermtime, sigchin);
	return (Cc + Cs*pi*pow(radtherm, 2)/sigchi2) * time;
}


/*============= MAIN  ================*/
int main(){
/*Cross-Sections and masses*/
	double DMmass = 0;
	double powsigchin =0;
	double powsigchi2 =0;
	cout << '\n';
	cout << "DM-Nucleon cross-section[order, cm^2]: "; cin >> powsigchin;
	cout << "DM-DM cross-section[order, cm^2]: "; cin >> powsigchi2;
	double SIGCHIN = pow(10., powsigchin);
	double SIGCHI2 = pow(10., powsigchi2);
	cout << "DM mass[GeV]: "; cin >> DMmass;
	double DMmass0 = DMmass;

	if (SIGCHI2/DMmass >= 2e-24){
		while(SIGCHI2/DMmass >=  2e-24){
			cout << '\n' << "Combination Breaks Bullet Cluster... Pick new ones" << '\n';
			cout << "DM-DM cross-section[order, cm^2]: "; cin >> powsigchi2;
			SIGCHI2 = pow(10., powsigchi2);
			cout << "DM mass[GeV]: "; cin >> DMmass;
		}
	}

	double ReducedMass = DMmass * Nmass/(DMmass + Nmass);
	double CaptureRate = 9.19e22 * (SIGCHIN/1e-55) / DMmass * rhoDM;		//Capture rate using these parameters
	double SelfCapture = 1.06e-3 * (SIGCHI2/1e-24) / DMmass * rhoDM;		//Self-capture rate using these parameters

/*Core Density */
	double nbee = 3*MNS/(Nmass*4*pi*pow(RNS,3));
	cout << "nB(pF): " << nbee << " ... nB:" << CoreDens << '\n';
	double TIME = 10;
	ofstream myfile;
	ofstream myfile2;
	ofstream myfile3;
	myfile.open("nchi.dat");
	myfile2.open("Eboundcheck.dat");
	myfile3.open("geotime.dat");
	
	double DMNO = 0;
	double DMNOplus = 0;
	
/* CHANDRASEKHAR */
	double chandpref = pow( (9./(32.*pi*pi)), 1./3.)*Planckbar*clight*5./3. /Gnewton;
	double Nchichand = pow( (chandpref/(DMmass*DMmass)), 3./2.);
	double Nchichandf = 1.8e51*pow((100/DMmass), 3.);
	double Nchichandb = 1.5e34*pow((100/DMmass), 2.);
/*TIMESCALES */
	double TGEO = tgeo(DMmass, ReducedMass, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
	double TTHERM = ttherm(DMmass, Nmass, ReducedMass, SIGCHIN);
	double times[4] = {TGEO, TTHERM, T1, TMAX};
	double t1, t2, t3 ,t4;
	int counter = 0;

	vector<double> timevec (times, times+4);
	double RTH = radius(DMmass, ReducedMass, TTHERM, SIGCHIN);
	double RTHEST = sqrt(9.* kBoltzmann * TNS/(4.*pi*Gnewton*CoreDens*DMmass));
	double RTHEST2 = 334*sqrt(Nmass/DMmass);
	cout << "containment time: " << T1 << '\n'
		<< "thermalization time: " << TTHERM << '\n'
		<< "goemetric time: " << TGEO << '\n'	
		<< "thermal. radius 1: " << RTH << '\n'	
		<< "thermal. radius 2: " << RTHEST << '\n'	
		<< "thermal. radius 3: " << RTHEST2 << '\n';	
	
/* RANK TIMES */
	stable_sort(timevec.begin(), timevec.end());
	t1 = timevec[0];
	t2 = timevec[1];
	t3 = timevec[2];
	t4 = timevec[3];
	double tprev = 0;
	while(TIME <= 1e11){
		if (TIME < t1){
			DMNOplus = 0;
		}
 
		else if (TIME >= t1 && TIME < t2){
			if(t1 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t1 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, Nmass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t1 == TMAX){
				DMNOplus = 0;
				cout << "DONE" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
			}
		}

		else if (TIME >= t2 && TIME < t3){
			if(t2 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t2 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, Nmass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t2 == TMAX){
				DMNOplus = 0;
				cout << "DONE" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
			}
		}


		else if (TIME >= t3 && TIME <t4){
			if(t3 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t3 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, Nmass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
			}
			else if(t3 == TMAX){
				DMNOplus = 0;
				cout << "DONE" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
			}
		}

		else if (TIME >= t4){
			DMNOplus = 0;
			if(counter == 0){
				cout << "DONE" << '\n';
				counter++;
			}
		}

	
		DMNO += DMNOplus;
		double Rchi = radius(DMmass, ReducedMass, TIME, SIGCHIN)*1e13;		//DM sphere radius in fm
		double Volumechi = 4./3. *pi * pow(Rchi, 3.);
		double DMDENS = DMNO/(Volumechi);			//DM number density in inverse fm^3
		myfile << TIME << " " << DMNO << " " << DMDENS << " " << SIGCHIN << " " << SIGCHI2 << " " << DMmass << " " << Nchichandf << " " << Nchichandb << '\n';
		tprev = TIME;
		TIME *= 1.1;
	}




/*================================================= OTHER TESTS =====================================================*/

	while(DMmass <= 1e4){
		double Einf = 2*Nmass*DMmass/pow((DMmass + Nmass), 2) * ( 1 - sqrt(1 - 2*Gnewton*MNS/(RNS*clight*clight)));
		double Emin = 1./3.*v2/(clight*clight);
		myfile2 << DMmass << " " << Emin << " " << Einf << '\n';
		DMmass *= 10;
	}

	TIME = 3.5; 
	ReducedMass = DMmass0 * Nmass/(DMmass0 + Nmass);
	while(TIME <= 1e10){
		double rad = radius(DMmass0, ReducedMass, TIME, SIGCHIN)/RNS;
		double rfid = sqrt(Nchi0(CaptureRate, SelfCapture, TIME)* SIGCHI2/(pi * RNS *RNS));

		myfile3 << TIME-3.5 << " " << rad << " " << rfid << '\n';	
//		cout << DMmass0 << ' ' << SIGCHIN << '\n';
		TIME *= 1.1;

	}
	myfile3.close();
	myfile2.close();
	myfile.close();
	return 0;
}

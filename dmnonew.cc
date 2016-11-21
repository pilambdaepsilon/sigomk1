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

/*============== UNIT CONVERSION natural units where c = hBAR = kB  = 1 ===========*/
double pi = M_PI;
const double convGeVtokg = 1.78e-27;
const double convGeVtoinvcm = 5.06e13;
const double convGeVtoinvs = 1.52e24;
const double convGeVtoK =1.16e13;

const double convcmtoinvGeV = 5.06e13;
const double convstoinvGeV = 1.52e24;

const double convkgtoGeV = 5.62e26;
const double convinvcmtoGeV = 1.98e-14;
const double convinvstoGeV = 6.58e-25;
const double convKtoGeV =8.62e-14;
const double convGeVtoinvFM = 0.0008065*(2*pi)*1e-3;			

/*===================== CONSTANTS ==========================*/
const double GNewton = 6.67408e-11;					// m^3/(kg*s^2)
const double kBoltzmann = 1.3806e-23;					// m^2*kg/(s^2*K)
const double clight = 3e8;						// m/s
const double Planckbar = (6.626e-34)/(2*pi);				// m^2*kg/s

const double GNewton0 = (6.67408e-11)*1e6*pow(convcmtoinvGeV,3.)/(convkgtoGeV*pow(convstoinvGeV,2.));			// natural, 1/M_planck^2
const double kBoltzmann0 = 1.3806e-23*1e4*pow(convcmtoinvGeV,2.)*convkgtoGeV/(pow(convstoinvGeV,2.)*convKtoGeV);	// natural, kB = 1
const double clight0 = 3e8*1e2*convcmtoinvGeV/convstoinvGeV;								// natural, c = 1
const double Planckbar0 = Planckbar*1e4*pow(convcmtoinvGeV,2)*convkgtoGeV/convstoinvGeV;				// natural, hBar = 1


/*====================== EMPERICAL NS PARAMETERS =======================*/
const double Msolar = 1.98892e30*convkgtoGeV;						// mass of Sol in GeV
double TNS = 1e5*convKtoGeV;								// NS temperature in K->GeV
double TDM = TNS;
double RNS = 10.6e5*convcmtoinvGeV;							// NS radius in cm->GeV^-1
double MNS = 1.44*Msolar;								// NS mass in GeV
double vbar = 220e5*convcmtoinvGeV/convstoinvGeV;					// mean thermal DM speed in cm/s->fraction of c - MB distrib
double v2 = vbar*vbar;									// the square pops up more frequently
const double Nmass = 0.94;								// Neutron mass in GeV
const double rhoDM = 1./pow(convcmtoinvGeV,3.);						// GeV/cm^3 -> GeV^4 - mean DM density at 1 kpc from Galactic center
const double rho_baryons = 5.7e11*convkgtoGeV/pow(convcmtoinvGeV, 3.);			// mean NS baryon mass density in kg/cm^3->GeV^4
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
	double CaptureRate = 9.19e22 * (SIGCHIN/1e-55) /DMmass;				//Capture rate using these parameters
	double SelfCapture = 1.06e-3 * (SIGCHI2/1e-24) /DMmass;				//Self-capture rate using these parameters

/*Core Density */
	double nbee = 3*MNS/(Nmass*4*pi*pow(RNS,3));
	cout << "nB(pF): " << nbee << " ... nB:" << CoreDens << '\n';
	double TIME = 1;
	ofstream myfile;
	ofstream myfile2;
	ofstream myfile3;
	myfile.open("nchi.dat");
	myfile2.open("Eboundcheck.dat");
	myfile3.open("geotime.dat");
	
	double DMNO = 0;
	double DMNOplus = 0;
	
/* CHANDRASEKHAR */
	double chandpref = pow( (9./(32.*pi*pi)), 1./3.)*Planckbar*clight*5./3. /GNewton;
	double Nchichand = pow( (chandpref/(DMmass*DMmass)), 3./2.);
	double Nchichandf = 1.8e51*pow((100/DMmass), 3.);
	double Nchichandb = 1.5e34*pow((100/DMmass), 2.);
/*TIMESCALES */
	double TCONTAIN = tcon(DMmass, SIGCHIN);
	double TGEO = tgeo(DMmass, ReducedMass, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
	double TTHERM = ttherm(DMmass, Nmass, ReducedMass, SIGCHIN);
	double times[4] = {TGEO, TTHERM, TCONTAIN, TMAX};
	double t1, t2, t3 ,t4;
	int counter = 0;

	vector<double> timevec (times, times+4);
	double RTH = radius(DMmass, ReducedMass, TTHERM, SIGCHIN)/convcmtoinvGeV;
	double RTHEST = sqrt(9.* TNS/(4.*pi*GNewton0*CoreDens*Nmass*DMmass))/convcmtoinvGeV;
	double RTHEST2 = 334*sqrt(Nmass/DMmass);
	cout << "containment time: " << TCONTAIN << '\n'
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
	cout << "t1: " << t1 << '\n'
		<< "t2: " << t2 << '\n'
		<< "t3: " << t3 << '\n'	
		<< "t4: " << t4 << '\n';	
	double tprev = 0;
	int count = 0;
	while(TIME <= TMAX){
		if (TIME < t1){
			DMNOplus = 0;
		}
 
		else if (TIME >= t1 && TIME < t2){
			if(t1 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 0){
				cout << "T1: Geometric" << '\n';
				count = 1;}
			}
			else if(t1 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 0){
				cout << "T1: Thermalization" << '\n';
				count = 1;}
			}
			else if(t1 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t1 = 10^10 years." << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 0){
				cout << "T1: Regular Evolution" << '\n';
				count = 1;}
			}
		}

		else if (TIME >= t2 && TIME < t3){
			if(t2 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 1){
				cout << "T2: Geometric" << '\n';
				count = 2;}
			}
			else if(t2 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 1){
				cout << "T2: Thermalization" << '\n';
				count = 2;}
			}
			else if(t2 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t2 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 1){
				cout << "T2: Regular Evolution" << '\n';
				count = 2;}
			}
		}


		else if (TIME >= t3 && TIME <t4){
			if(t3 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 2){
				cout << "T3: Geometric" << '\n';
				count = 3;}
			}
			else if(t3 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 2){
				cout << "T3: Thermalization" << '\n';
				count = 3;}
			}
			else if(t3 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t3 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 2){
				cout << "T3: Regular Evolution" << '\n';
				count = 3;}
			}
		}

		else if (TIME >= t4){
			DMNOplus = 0;
			if(counter == 0){
				cout << "DONE, t4 = 10^10 years" << '\n';
				counter++;
			}
		}

	
		DMNO += DMNOplus;
		double Rchi = radius(DMmass, ReducedMass, TIME, SIGCHIN)*1e13;		//DM sphere radius in fm
		double Volumechi = 4./3. *pi * pow(Rchi, 3.);
		double DMDENS = DMNO/(Volumechi);			//DM number density in inverse fm^3
		myfile << TIME << " " << DMNO << " " << DMDENS << " " << SIGCHIN << " " << SIGCHI2 << " " << DMmass << " " << Nchichandf << " " << Nchichandb << '\n';
		tprev = TIME;
		TIME *= 1.01;
	}




/*================================================= OTHER TESTS =====================================================*/

	while(DMmass <= 1e4){
		double Einf = 2*Nmass*DMmass/pow((DMmass + Nmass), 2) * ( 1 - sqrt(1 - 2*GNewton*MNS/(RNS*clight*clight)));
		double Emin = 1./3.*v2/(clight*clight);
		myfile2 << DMmass << " " << Emin << " " << Einf << '\n';
		DMmass *= 10;
	}

	TIME = 3.5; 
	ReducedMass = DMmass0 * Nmass/(DMmass0 + Nmass);
	while(TIME <= 1e10){
		double rad = radius(DMmass0, ReducedMass, TIME, SIGCHIN)/RNS;
		double rfid = sqrt(Nchi0(CaptureRate, SelfCapture, TIME)* SIGCHI2*pow(convcmtoinvGeV,2.)/(pi * RNS *RNS));

		myfile3 << TIME-3.5 << " " << rad << " " << rfid << '\n';	
//		cout << DMmass0 << ' ' << SIGCHIN << '\n';
		TIME *= 1.001;

	}
	myfile3.close();
	myfile2.close();
	myfile.close();
	return 0;
}

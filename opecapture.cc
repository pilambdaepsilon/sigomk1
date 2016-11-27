#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"


double coeff(int counter, double mDM, double sigchi2, double w, double vT){
	sigchi2 *= pow(convcmtoinvGeV, 2.);			//use natural units for the cross-section
	if(counter == 1){
		return sqrt(sigchi2 * 16. * pi/pow(mDM, 1.));
	}
	else if(counter == 2){
		return 0.0;
	}
	else if(counter == 3){
		return sqrt(sigchi2 * 128. * pi/(mDM*pow(vT*w, 2.)));
	}
	else if(counter == 4){
		return sqrt(sigchi2 * 64. * pi/(mDM));
	}
	else if(counter == 5){
		return sqrt(sigchi2 * 96. * pi/(mDM*pow(vT*w, 2.)));
	}
	else if(counter == 6){
		return sqrt(sigchi2 * 576. * pi/(mDM*pow(w*w, 2.)));
	}
	else if(counter == 7){
		return sqrt(sigchi2 * 64. * pi/(mDM*pow(vT, 2.)));
	}
	else if(counter == 8){
		return sqrt(sigchi2 * 48. * pi/(mDM*pow(vT, 2.)));
	}
	else if(counter == 9){
		return sqrt(sigchi2 * 96. * pi/(mDM*pow(vT, 2.)));
	}
	else if(counter == 10){
		return sqrt(sigchi2 * 192. * pi/(mDM*pow(w, 2.)));
	}
	else if(counter == 11){
		return sqrt(sigchi2 * 128.  * pi/(mDM*pow(w, 2.)));
	}
	else if(counter == 12){
		return sqrt(sigchi2 * 96. * pi/(mDM*pow(w, 2.)));
	}
	else if(counter == 13){
		return sqrt(sigchi2 * 384. * pi/(mDM*pow(vT*w, 2.)));
	}
	else if(counter == 14){
		return sqrt(sigchi2 * 384. * pi/(mDM*pow(vT*w, 2.)));
	}
	else if(counter == 15){
		return sqrt(sigchi2 * 576. * pi/(mDM*pow(vT*w*w, 2.)));
	}
	else{
		return 0.0;
	}
}
double crosssection(int counter, double mDM, double coeff, double w, double vT){
	if(counter == 1){
		return pow(mDM, 2.)/(64.*pi) * 4.*coeff;
	}
	else if(counter == 2){
		return 0;
	}
	else if(counter == 3){
		return pow(mDM*w*vT, 2.)/(384.*pi) * 3.*coeff;
	}
	else if(counter == 4){
		return pow(mDM, 2.)/(64.*pi) * coeff;
	}
	else if(counter == 5){
		return pow(mDM*w*vT, 2.)/(384.*pi) * 4.*coeff;
	}
	else if(counter == 6){
		return pow(mDM*w*w, 2.)/(576.*pi) *coeff;
	}
	else if(counter == 7){
		return pow(mDM*vT, 2.)/(192.*pi) * 3.*coeff;
	}
	else if(counter == 8){
		return pow(mDM*vT, 2.)/(192.*pi) * 4.*coeff;
	}
	else if(counter == 9){
		return pow(mDM*vT, 2.)/(192.*pi) * 2.*coeff;
	}
	else if(counter == 10){
		return pow(mDM*w, 2.)/(384.*pi) * 2.*coeff;
	}
	else if(counter == 11){
		return pow(mDM*w, 2.)/(384.*pi) * 3.*coeff;
	}
	else if(counter == 12){
		return pow(mDM*w, 2.)/(384.*pi) * 4.*coeff;
	}
	else if(counter == 13){
		return pow(mDM*w*vT, 2.)/(384.*pi) * coeff;
	}
	else if(counter == 14){
		return pow(mDM*w*vT, 2.)/(384.*pi) * coeff;
	}
	else if(counter == 15){
		return pow(mDM*w*w*vT, 2.)/(576.*pi) *coeff;
	}
	else{
		return 0.0;
	}
}
double capture(double mDM, double rx, double sigchi2){
	double vd = 220e5*convcmtoinvGeV/convstoinvGeV;
	double vo = 4e-3;
	return rx/mDM * pow(vo,2.)/ vd * sigchi2;
}

int main(){
	double DMmass = 10;		//GeV
	double Nmass = 1;
	double ReducedMass = DMmass*Nmass/(DMmass+Nmass);
	double ReducedMassDM = DMmass*DMmass/(DMmass+DMmass);
	double DMdens = 0.4/pow(convcmtoinvGeV,3.);		//GeV/cm^3 ->GeV^4
	double SIGCHI2 = 1.78e-25;	//cm^2 ->>>>>>> constraint from N-body simulations on DM self-interaction x-sec
//	double SIGCHI2 = 2.23e-24;	//cm^2 ->>>>>>> constraint from Bullet Cluster on DM self-interaction x-sec
	double W = 1e8;		//cm/s
	double VT = W;		//cm/s
//	double VT = W*(1 - DMmass*DMmass/(4*ReducedMass*ReducedMass));		//cm/s
	double dummyc = 0;
	double Win = 0;
	cout << "vT: "; cin >> Win; cout <<'\n';
	double VTin = Win;		//cm/s
//	double VTin = Win*(1 - DMmass*DMmass/(4*ReducedMass*ReducedMass));		//cm/s
	cout << '\n' << "w: " << W << " vT: " << VT << '\n';
	ofstream myfile;
	myfile.open("coeffs.dat");

	ofstream myfile2;
	myfile2.open("xsecs.dat");

	ofstream myfile3;
	myfile3.open("rates.dat");

	VT *= convcmtoinvGeV/convstoinvGeV;	// as a fraction of c
	W *= convcmtoinvGeV/convstoinvGeV;	// as a fraction of c
	Win *= convcmtoinvGeV/convstoinvGeV;	// as a fraction of c
	VTin *= convcmtoinvGeV/convstoinvGeV;	
	double Coefficients[16] = {0.0};
	double Xsections[16] = {0.0};
	double Capturerates[16] = {0.0};

	while(DMmass <= 2e3){
		/* ================= GET COEFFICIENTS =============== */
		for(int i = 0; i<16; i++){
			Coefficients[i] = coeff(i+1, DMmass, SIGCHI2, W, VT);
		}
		myfile << DMmass << ' ';
		for(int i = 0; i<15; i++){
			myfile << Coefficients[i] << ' ';
		}
			myfile << '\n';

		/* ================ GET CROSS-SECTIONS ============== */
		for(int i = 0; i<16; i++){
			Xsections[i] = crosssection(i+1, DMmass, Coefficients[i]*Coefficients[i], Win, VTin)*pow(convinvGeVtocm,2.);
		}
		myfile2 << DMmass << ' ';
		for(int i = 0; i<15; i++){
			myfile2 << Xsections[i] << ' ';
		}
			myfile2 << 2.23e-24*DMmass <<  '\n';
		if (DMmass == 10){
			cout << '\n' << '\n';
				cout << "=============== X-SEC @ 10 GeV mass ===============" << '\n';
				cout << DMmass << ' ';
				for(int i = 0; i<15; i++){
					cout << Xsections[i] << ' ';
				}
				cout << '\n' << "==================================================";
			cout << '\n' << '\n';
		}
			
		/* ================ GET CAPTURE RATE ================ */
		for(int i = 0; i<16; i++){
			Capturerates[i] = capture(DMmass, DMdens, Xsections[i]*pow(convcmtoinvGeV,2.))*convGeVtoinvs;
		}
		myfile3 << DMmass << ' ';
		for(int i = 0; i<15; i++){
			myfile3 << Capturerates[i] << ' ';
		}
			myfile3 << '\n';

		/*
		myfile << DMmass << ' ' <<
		Xsections[0] + Xsections[3] << ' ' <<
		Xsections[6] + Xsections[7]  + Xsections[8] + Xsections[9] + Xsections[10] + Xsections[11]  << ' ' <<
		Xsections[2] + Xsections[4] + Xsections[5] + Xsections[12] + Xsections[13] << ' ' <<
		Xsections[14] << '\n';
		*/

		DMmass *= 1.1;
	}
	myfile.close();		
	myfile2.close();	
	myfile3.close();	

	cout << DMmass << ' ' <<
	Xsections[0] + Xsections[3] << ' ' <<
	Xsections[6] + Xsections[7]  + Xsections[8] + Xsections[9] + Xsections[10] + Xsections[11]  << ' ' <<
	Xsections[2] + Xsections[4] + Xsections[5] + Xsections[12] + Xsections[13] << ' ' <<
	Xsections[14] << '\n' << '\n';

		for(int i = 0; i<15; i++){
			cout << Coefficients[i] << ' ';
		}
	cout << '\n' << '\n';
//		for(int i = 0; i<15; i++){
//			cout << Xsections[i] << ' ';
//		}
//	cout << '\n' << '\n';
		for(int i = 0; i<15; i++){
			cout << Capturerates[i] << ' ';
		}
	cout << '\n' << '\n';

	return 0;
}

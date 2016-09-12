#include "sigofunctions.hh"

double Int_rhoS(double kmin, double kmax, double MEFF, int calls){
	gsl_function RhoScalar;
	RhoScalar.function = &rhoS;
	RhoScalar.params = &MEFF;

	gsl_integration_workspace *WR = gsl_integration_workspace_alloc(calls);
	double relerr=1e-7;

	double RHOSC = 0;
	double RHOERR = 0;
	
	gsl_integration_qags(&RhoScalar, kmin, kmax, 0., relerr, calls, WR, &RHOSC, &RHOERR);
	gsl_integration_workspace_free(WR);

	return RHOSC;
} 

double Int_Ekinetic(double kmin, double kmax, double mst, int calls){
	gsl_function EK;

	EK.function = &EKIN;
	EK.params = &mst;

	gsl_integration_workspace * WE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double EKin = 0;
	double EKinERR = 0;

	gsl_integration_qags(&EK, kmin, kmax, 0., relerr0, calls, WE, &EKin, &EKinERR);
	gsl_integration_workspace_free(WE);
	
	return EKin;
}
double Int_eye1(double kmin, double kmax, int calls){
	gsl_function EYE1;

	EYE1.function = &I1;
	EYE1.params = &mstar;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Ione = 0;
	double Ioneerr = 0;

	gsl_integration_qags(&EYE1, kmin, kmax, 0., relerr0, calls, WEYE, &Ione, &Ioneerr);
	gsl_integration_workspace_free(WEYE);
	
	return Ione;
}

double Int_eye2(double kmin, double kmax, int calls){
	gsl_function EYE2;

	EYE2.function = &I2;
	EYE2.params = &mstar;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Itwo = 0;
	double Itwoerr = 0;

	gsl_integration_qags(&EYE2, kmin, kmax, 0., relerr0, calls, WEYE, &Itwo, &Itwoerr);
	gsl_integration_workspace_free(WEYE);
	
	return Itwo;
}

double Int_eye3(double kmin, double kmax, int calls){
	gsl_function EYE3;

	EYE3.function = &I3;
	EYE3.params = &mstar;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Ithree = 0;
	double Ithreeerr = 0;

	gsl_integration_qags(&EYE3, kmin, kmax, 0., relerr0, calls, WEYE, &Ithree, &Ithreeerr);
	gsl_integration_workspace_free(WEYE);
	
	return Ithree;
}


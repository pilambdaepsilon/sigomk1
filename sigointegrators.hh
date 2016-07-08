#include "sigofunctions.hh"

double Int_rhoS(double kmin, double kmax, int calls){
	gsl_function RhoScalar;
	struct arg_params Rhoparams = {mass};
	RhoScalar.function = &rhoS;
	RhoScalar.params = &Rhoparams;

	gsl_integration_workspace *WR = gsl_integration_workspace_alloc(calls);
	double relerr=1e-7;

	double RHOSC = 0;
	double RHOERR = 0;
	
	gsl_integration_qags(&RhoScalar, kmin, kmax, 0., relerr, calls, WR, &RHOSC, &RHOERR);
	gsl_integration_workspace_free(WR);

	return RHOSC;
} 

double Int_Ekinetic(double kmin, double kmax, int calls){
	gsl_function EK;
	struct arg_params EParams = {mass};

	EK.function = &goko;
	EK.params = &EParams;

	gsl_integration_workspace * WE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double EKin = 0;
	double EKinERR = 0;

	gsl_integration_qags(&EK, kmin, kmax, 0., relerr0, calls, WE, &EKin, &EKinERR);
	gsl_integration_workspace_free(WE);
	
	return EKin;
}


#include "SigOmRhoInt.hh"
double Pi = pi;
double Intecheck1 = 0;
double Intecheck2 = 0;
double Intecheck3 = 0;
double mstarcheck = 0;
struct Rparams{
	double Msig; double Mw;	double Mr; double Mn; double Me; double gs; double gw; double gr;
	double bethe; double cethe; double I3; double Jn; double Jp; double Je;	double ban; double bap;		
	double Qe; double Qp; double nB; double Easymm;	double Kcompress; double BperA;	double gamma; double asymmetryfactor;
	double Integ1; double Integ2; double Integ3;
};

//DEFINE FUNCTION FOR SATURATION CONSTANTS
int GM1_f(const gsl_vector *x, void *params, gsl_vector *f){
	//SET PARAMETERS
	double Mn = ((struct Rparams *) params)->Mn;
	double Me = ((struct Rparams *) params)->Me;
	double I3 = ((struct Rparams *) params)->I3;
	double Jn = ((struct Rparams *) params)->Jn;
	double Jp = ((struct Rparams *) params)->Jp;
	double Je = ((struct Rparams *) params)->Je;
	double Qe = ((struct Rparams *) params)->Qe;
	double Qp = ((struct Rparams *) params)->Qp;
	double nB = ((struct Rparams *) params)->nB;
	double Easymm = ((struct Rparams *) params)->Easymm;
	double Kcompress = ((struct Rparams *) params)->Kcompress;
	double BperA = ((struct Rparams *) params)->BperA;
	double gamma = ((struct Rparams *) params)->gamma;
	double asymmetryfactor = ((struct Rparams *) params)->asymmetryfactor;
	
	//SET VARIABLES TO BE SOLVED FOR
	double gsigma = gsl_vector_get(x,0);
	double gomega = gsl_vector_get(x,1);
	double grho = gsl_vector_get(x,2);
	double asigma = gsl_vector_get(x,3);
	double aomega = gsl_vector_get(x,4);
	double arho = gsl_vector_get(x,5);
	double b = gsl_vector_get(x,6);
	double c = gsl_vector_get(x,7);

	//SET FUNCTIONS TO BE SOLVED
	double mstar = Mn*gamma;
	double kF = pow((3*Pi*Pi*nB*0.5),1./3); 
	double kFactor = sqrt(kF*kF + mstar*mstar);	
	int calls = 1000;
//	double Integral1 = Int_eye1(0, kF, mstar, calls); 
	double Integral1 = eye1(kF, mstar); 
//	double Integral2 = Int_eye2(0, kF, mstar, calls); 
	double Integral2 = eye2(kF, mstar); 
//	double Integral3 = Int_eye3(0, kF, mstar, calls); 
	double Integral3 = eye3(kF, mstar); 
	double Pref7 = 1./asigma + 2*b*Mn*gsigma + 3*c*pow(gsigma,2.) + 2/(Pi*Pi)*Integral3;
	
	double f0 = gomega - aomega*nB;												//*
	double f1 = grho + 0.5*arho*asymmetryfactor*nB;										//*
	double f2 = gsigma + b* Mn*asigma*pow(gsigma, 2.) + c*asigma*pow(gsigma,3.) - 2./(Pi*Pi)*asigma*Integral1;		//*
	double f3 = nB*aomega - Mn - BperA + kFactor;
	double f4 = gsigma - (1 - gamma)*Mn;
 	double f5 = nB*(Mn +BperA) - 2./(Pi*Pi)*Integral2 - 1./3*b*Mn*pow(gsigma,3.) - 0.25*c*pow(gsigma,4.) - 0.5*gsigma*gsigma/asigma - 0.5*gomega*gomega/aomega;
	double f6 = Easymm - arho*nB/(8.) - kF*kF/(6*kFactor);
 	double f7 = Kcompress - 6*aomega/(Pi*Pi)*pow(kF,3.) - 3*kF*kF/kFactor + 6*pow(kF,3.)/(Pi*Pi)*mstar*mstar/pow(kFactor, 2.)/Pref7; 
	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	gsl_vector_set(f, 2, f2);
	gsl_vector_set(f, 3, f3);
	gsl_vector_set(f, 4, f4);
	gsl_vector_set(f, 5, f5);
	gsl_vector_set(f, 6, f6);
	gsl_vector_set(f, 7, f7);
	
	return GSL_SUCCESS;
}

//DEFINE FUNCTION FOR EVOLUTION with nB
int GM2_f(const gsl_vector *x, void *params, gsl_vector *f){
	//SET PARAMETERS
	double Mn = ((struct Rparams *) params)->Mn;
	double Me = ((struct Rparams *) params)->Me;
	double as = ((struct Rparams *) params)->gs;
	double aw = ((struct Rparams *) params)->gw;
	double ar = ((struct Rparams *) params)->gr;
	double bethe = ((struct Rparams *) params)->bethe;
	double cethe = ((struct Rparams *) params)->cethe;
	double I3 = ((struct Rparams *) params)->I3;
	double Jn = ((struct Rparams *) params)->Jn;
	double Jp = ((struct Rparams *) params)->Jp;
	double Je = ((struct Rparams *) params)->Je;
	double Qe = ((struct Rparams *) params)->Qe;
	double Qp = ((struct Rparams *) params)->Qp;
	double nB = ((struct Rparams *) params)->nB;
	double asymmetryfactor = ((struct Rparams *) params)->asymmetryfactor;

	//SET VARIABLES TO BE SOLVED FOR
	double gsigma = gsl_vector_get(x,0);
	double gomega = gsl_vector_get(x,1);
	double grho = gsl_vector_get(x,2);
	double mun = gsl_vector_get(x,3);
	double mup = gsl_vector_get(x,4);
	double mue = gsl_vector_get(x,5);

	//SET FUNCTIONS TO BE SOLVED
	double mstar = Mn - gsigma;
	double kF = pow((3*Pi*Pi*nB*0.5),1./3); 
	double kFactor = sqrt(kF*kF + mstar*mstar);	
	double kFasymn = kF*pow((1 + asymmetryfactor)*0.5,1./3);
	double kFasymp = kF*pow((1 - asymmetryfactor)*0.5,1./3);
	double Integral1 = eye1(kF, mstar); 
	double Integral2 = eye2(kF, mstar); 
	double Integral3 = eye3(kF, mstar); 
	
	double f0 = gomega - aw*nB;
	double f1 = grho + 0.5*ar*nB*asymmetryfactor;
	double f2 = gsigma + bethe* Mn*as*pow(gsigma, 2.) + cethe*as*pow(gsigma,3.) - 2./(Pi*Pi)*as*Integral1;
	double f3 = mun - sqrt(kFasymn*kFasymn + mstar*mstar) - gomega + 0.5*grho;
	double f4 = mup - sqrt(kFasymp*kFasymp + mstar*mstar) - gomega - 0.5*grho;
	double f5 = mue - mun + mup;
	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	gsl_vector_set(f, 2, f2);
	gsl_vector_set(f, 3, f3);
	gsl_vector_set(f, 4, f4);
	gsl_vector_set(f, 5, f5);
	
	return GSL_SUCCESS;
}


int print_state(size_t iter, gsl_multiroot_fsolver *S){
cout << "iter = " << iter << '\n' << 
	" x = " << "gsigma: " << gsl_vector_get(S->x, 0)/MeVtoinvFM << ", gomega: " << gsl_vector_get(S->x, 1)/MeVtoinvFM << ", grho: " << gsl_vector_get(S->x, 2)/MeVtoinvFM 
	<< ", asigma: " << gsl_vector_get(S->x, 3) << ", aomega: " <<
	gsl_vector_get(S->x, 4) << ", arho: " << gsl_vector_get(S->x, 5) << ", b: " << gsl_vector_get(S->x, 6) << ", c: " << gsl_vector_get(S->x, 7) << '\n' << 
	" F(x) = " << gsl_vector_get(S->f, 0) << ", " << gsl_vector_get(S->f, 1) << ", " << gsl_vector_get(S->f, 2) << ", " << gsl_vector_get(S->f, 3) << ", " <<
	gsl_vector_get(S->f, 4) << ", " << gsl_vector_get(S->f, 5) << ", " << gsl_vector_get(S->f, 6) << ", " << gsl_vector_get(S->f, 7) << '\n' <<
	 "============================================================" << '\n'; 
}
int print_state2(size_t iter, gsl_multiroot_fsolver *S2){
cout << "iter = " << iter << '\n' << 
	" x = " << "gsigma: " << gsl_vector_get(S2->x, 0)/MeVtoinvFM << ", gomega: " << gsl_vector_get(S2->x, 1)/MeVtoinvFM << ", grho: " <<
       	gsl_vector_get(S2->x, 2)/MeVtoinvFM << ", MuN: " << gsl_vector_get(S2->x, 3)/MeVtoinvFM << ", MuP: " << gsl_vector_get(S2->x, 4)/MeVtoinvFM << ", MuE: " << gsl_vector_get(S2->x, 5)/MeVtoinvFM << '\n' <<
	" F(x) = " << gsl_vector_get(S2->f, 0) << ", " << gsl_vector_get(S2->f, 1) << ", " << gsl_vector_get(S2->f, 2) << ", "  << gsl_vector_get(S2->f,3) << ", " <<
	gsl_vector_get(S2->f,4) << ", " << gsl_vector_get(S2->f, 5) << '\n' <<
	 "============================================================" << '\n'; 
}


int main(){
	double E0 = 0;
	int calls = 1000;
	double BaryonDensity = 0.153;

/*================================== DEFINE PARAMETERS TO PASS TO SOLVER ==============================*/
	double ScalarMass = 550*MeVtoinvFM;
	double VectorMass = 783*MeVtoinvFM;
	double IsoVectorMass = 800*MeVtoinvFM;
	double NMass = 938*MeVtoinvFM;
	double ElectronMass = 0.511*MeVtoinvFM;
	double ScalarField = 0.0;
	double VectorField = 0.0;
	double IsoVectorField = 0.0;
	double B = 0.0;
	double C = 0.0;
	double IsospinN = 0.5;
	double SpinN = 0.5;
	double SpinP = 0.5;
	double SpinE = 0.5;
	double BaryonNumberN = 1;
	double BaryonNumberP = 1;
	double ChargeP = 1;
	double ChargeE = -1;
	double AsymmetryEnergy = 32.5*MeVtoinvFM;
	double Compressibility = 240*MeVtoinvFM;
	double BindingPerNucleon0 = -16.3*MeVtoinvFM;
	double MstarPerMn0 = 0.78;
	double AlphaAsymmetry = 0.;			// Aplha = 0 is for symmetric matter 
	double I1 = 0;
	double I2 = 0;
	double I3 = 0;
	
	double kF = pow((3*Pi*Pi*BaryonDensity*0.5),1./3); 
	I1 = Int_eye1(0, kF, MstarPerMn0, calls); 
	I2 = Int_eye2(0, kF, MstarPerMn0, calls); 
	I3 = Int_eye3(0, kF, MstarPerMn0, calls); 

	struct Rparams p = {ScalarMass, VectorMass,  IsoVectorMass, NMass, ElectronMass, ScalarField,  VectorField, IsoVectorField, B, C, IsospinN, SpinN, SpinP, SpinE, BaryonNumberN,
		BaryonNumberP, ChargeE, ChargeP, BaryonDensity, AsymmetryEnergy, Compressibility, BindingPerNucleon0, MstarPerMn0, AlphaAsymmetry, I1, I2, I3};		//INITIALIZE PARAMETERS
/*=============================================================================================================*/

	ofstream DATER;
	DATER.open("rootE.dat");		//OPEN DATA FILES TO WRITE TO

/* ======================================= SET UP SOLVER & SOLVE SYSTEM AT SATURATION ================================================*/
	const gsl_multiroot_fsolver_type *T;	//DECLARE VARIABLE FOR TYPE OF SOLVER
	gsl_multiroot_fsolver *S;		//DECLARE NAME FOR SOLVER WORKSPACE

	int status;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
	size_t i, iter =0;

	const size_t n = 8;			//DIMENSIONALITY OF PROBLEM
				//  Msig,  Mw,   Mr,  Mn,   Me,  gs, gw, gr,   b,    c,    I3,  Jn,  Jp,  Je,  ban, bap, Qe,  Qp, nB             Easym, Kcompress   Eo 

	gsl_multiroot_function func = {&GM1_f, n, &p};	//MAKE A FUNCTION FOR SOLVER OUT OF ROSENBROCK, OF DIMENSION n, WITH PARAMS p
	double x_initial[8] = {1, 1, 1, 1, 1, 1, 1e-5, -1e-5};	//MAKE INITIAL GUESS FOR VARIABLES TO BE SOLVED FOR
	gsl_vector *x = gsl_vector_alloc(n);	//MAKE A VECTOR FOR SOLVER TO STORE VARIABLES
	
	gsl_vector_set(x, 0, x_initial[0]);
	gsl_vector_set(x, 1, x_initial[1]);	//SET VARIABLES TO BE SOLVED FOR TO INITIAL VALUES
	gsl_vector_set(x, 2, x_initial[2]);
	gsl_vector_set(x, 3, x_initial[3]);
	gsl_vector_set(x, 4, x_initial[4]);
	gsl_vector_set(x, 5, x_initial[5]);
	gsl_vector_set(x, 6, x_initial[6]);
	gsl_vector_set(x, 7, x_initial[7]);
	
	T = gsl_multiroot_fsolver_hybrids;	//MAKE THE TYPE OF SOLVER A HYBRID S SOLVER
	S = gsl_multiroot_fsolver_alloc(T, 8);	//MAKE THE WORKSPACE FOR SOLVER TYPE T OF DIMENSION 2
	gsl_multiroot_fsolver_set(S, &func, x);	//USE THE SOLVER TO SOLVE FUNCTIONS func AND FOR VARIABLES x (both have to be gsl_multiroots
	
		//print_state(iter, S);			//PRINT THE ORIGINAL STATE
	while(status = GSL_CONTINUE && iter < 1000){
		iter ++;
		status = gsl_multiroot_fsolver_iterate (S);	//ITERATE THE WORKSPACE UNTIL THE THINGS IS SOLVED

		if(status) break;				//CHECK IF IT'S WORKING. IF NOT, STOP
		status = gsl_multiroot_test_residual(S->f, 1e-7);	//SET TOLERANCE FOR CONVERGENCE
	}

/*========================== PULL CONSTANTS FROM SOLUTION =================================*/
	ScalarField = gsl_vector_get(S->x, 0); 
	VectorField = gsl_vector_get(S->x, 1);
	IsoVectorField = gsl_vector_get(S->x, 2);
	double ChemicalPotentialNeutron = 1.;
	double ChemicalPotentialProton = 1.;
	double ChemicalPotentialElectron = 2.;
	double Asigma = gsl_vector_get(S->x, 3);
	double Aomega = gsl_vector_get(S->x, 4);
	double Arho = gsl_vector_get(S->x, 5);
	B = gsl_vector_get(S->x, 6);
	C = gsl_vector_get(S->x, 7); 

/*================= SOLVE FOR FIELD AMPLITUDES WITH KNOWN CONSTANTS FROM SATURATION ========*/
	ScalarField = 0.1;
	VectorField = 1.;
	IsoVectorField = 0.5;
	BaryonDensity = 0.001;	
	while (BaryonDensity <= 2.0){
		AlphaAsymmetry = 0.178*log(40./(BaryonDensity+0.153));		//This is an approximate fit of the asymmetryfactor as a function of total Baryon Density
//		AlphaAsymmetry = 0.;
 		double NeutronDensity = 0.5*BaryonDensity*(1 + AlphaAsymmetry);
		double ProtonDensity = (BaryonDensity - NeutronDensity);
		
		const gsl_multiroot_fsolver_type *T2;		//DECLARE VARIABLE FOR TYPE OF SOLVER
		gsl_multiroot_fsolver *S2;		//DECLARE NAME FOR SOLVER WORKSPACE

		int status2;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
		size_t i2, iter2 =0;
	
		const size_t n2 = 6;			//DIMENSIONALITY OF PROBLEM
				//  Msig,  Mw,   Mr,  Mn,   Me,  gs, gw, gr,   b,    c,    I3,  Jn,  Jp,  Je,  ban, bap, Qe,  Qp, nB             Easym, Kcompress   Eo 
		struct Rparams p = {ScalarMass, VectorMass,  IsoVectorMass, NMass, ElectronMass, Asigma, Aomega, Arho, B, C, IsospinN, SpinN, SpinP, SpinE,
		       	BaryonNumberN, BaryonNumberP, ChargeE, ChargeP, BaryonDensity, AsymmetryEnergy, Compressibility, BindingPerNucleon0, MstarPerMn0, AlphaAsymmetry, I1, I2, I3};

		gsl_multiroot_function func2 = {&GM2_f, n2, &p};	//MAKE A FUNCTION FOR SOLVER OUT OF ROSENBROCK, OF DIMENSION n, WITH PARAMS p

		double y_initial[6] = {ScalarField, VectorField, IsoVectorField, 0.2, 0.1, 0.1};	

		gsl_vector *y = gsl_vector_alloc(n2);	//MAKE A VECTOR FOR SOLVER TO STORE VARIABLES
	
		gsl_vector_set(y, 0, y_initial[0]);
		gsl_vector_set(y, 1, y_initial[1]);	//SET VARIABLES TO BE SOLVED FOR TO INITIAL VALUES
		gsl_vector_set(y, 2, y_initial[2]);
		gsl_vector_set(y, 3, y_initial[3]);
		gsl_vector_set(y, 4, y_initial[4]);
		gsl_vector_set(y, 5, y_initial[5]);
	
		T2 = gsl_multiroot_fsolver_hybrids;	//MAKE THE TYPE OF SOLVER A HYBRID S SOLVER
		S2 = gsl_multiroot_fsolver_alloc(T, n2);	//MAKE THE WORKSPACE FOR SOLVER TYPE T OF DIMENSION n2
		gsl_multiroot_fsolver_set(S2, &func2, y);	//USE THE SOLVER TO SOLVE FUNCTIONS func AND FOR VARIABLES x (both have to be gsl_multiroots
	
		while(status2 = GSL_CONTINUE && iter2 < 1000){
			iter2 ++;
			status2 = gsl_multiroot_fsolver_iterate (S2);	//ITERATE THE WORKSPACE UNTIL THE THINGS IS SOLVED

			if(status2) break;				//CHECK IF IT'S WORKING. IF NOT, STOP
			status2 = gsl_multiroot_test_residual(S2->f, 1e-7);	//SET TOLERANCE FOR CONVERGENCE
		}

		ScalarField = gsl_vector_get(S2->x,0);
		VectorField = gsl_vector_get(S2->x,1);
		IsoVectorField = gsl_vector_get(S2->x,2);

		ChemicalPotentialNeutron = gsl_vector_get(S2->x,3);
		ChemicalPotentialProton = gsl_vector_get(S2->x,4);
		ChemicalPotentialElectron = gsl_vector_get(S2->x,5);

		double kFermiNeutron = pow((3*Pi*Pi*NeutronDensity*0.5),1./3); 
		double kFermiProton = pow((3*Pi*Pi*ProtonDensity*0.5),1./3); 
		double kFermiElectron = sqrt(ChemicalPotentialElectron*ChemicalPotentialElectron - ElectronMass*ElectronMass);
		double kFermi = pow((3*Pi*Pi*BaryonDensity*0.5),1./3); 
		double Mstar = NMass - ScalarField;
		double ScalarDensity = 2./(Pi*Pi)*(Int_eye1(0, kFermiNeutron, Mstar, calls) + Int_eye1(0, kFermiProton, Mstar, calls));
		double ScalarDensity2 = 2./(Pi*Pi)*Int_eye1(0, kFermi, Mstar, calls);


		double EOmega = 0.5*pow(VectorField, 2.)/Aomega;
		double EOmega2 = 0.5*Aomega*pow(BaryonDensity, 2.);

		double ERho = 0.5*IsoVectorField*IsoVectorField/Arho;
		double ERho2 = 1./8*Arho*pow((ProtonDensity - NeutronDensity), 2.);
	
		double ESigma = 0.5*pow(ScalarField,2.)/Asigma;
		double ESigma2 = 0.5*Asigma*pow(ScalarDensity2, 2.);
		double USigma =( (1./3.)*B*NMass*pow(ScalarField, 3.) + (0.25) * C*pow(ScalarField, 4.) ); 
 
		double EkinProton = 2./(Pi*Pi) * Int_eye2(0, kFermiProton, Mstar, calls);
		double EkinNeutron = 2./(Pi*Pi) * Int_eye2(0, kFermiNeutron, Mstar, calls);
		double Ekin = EkinNeutron + EkinProton;
		double Ekin2 = 2./(Pi*Pi)*Int_eye2(0, kFermi, Mstar, calls);
		double EElectron = 0.;
		double PElectron = 0.;
		if (ChemicalPotentialElectron >= ElectronMass){
			EElectron = 2./(Pi*Pi)*Int_eye2(0, kFermiElectron, ElectronMass, calls);
			PElectron = 1./(3*Pi*Pi)*ChemicalPotentialElectron*pow((ChemicalPotentialElectron*ChemicalPotentialElectron - ElectronMass*ElectronMass),3./2)
				- EElectron;
		}
		else{
		EElectron = 0.0;
		PElectron = 0.0;}

		double EnergyDensity = 
		 Ekin2
		
		+ EOmega
		+ ERho
		+ ESigma
		+ USigma;

		double BindingPerNucleon = (EnergyDensity/BaryonDensity - NMass);

		double Pressure = 
		EOmega	       
		+ ERho
		- ESigma
		- USigma
		+ 1./(3*Pi*Pi)*Int_eyeP(0, kFermiProton, Mstar, calls)  + 1./(3*Pi*Pi)*Int_eyeP(0, kFermiNeutron, Mstar, calls)
		;  

		DATER << BaryonDensity << " " << (BindingPerNucleon)/MeVtoinvFM << " " << Pressure/MeVtoinvFM << " " << ScalarDensity << " " << AlphaAsymmetry << " " << 
			Mstar/MeVtoinvFM << " " << VectorField/MeVtoinvFM << " " << -IsoVectorField/MeVtoinvFM << " " << ScalarField/MeVtoinvFM << " " << 
			ChemicalPotentialElectron/MeVtoinvFM << " " << (ChemicalPotentialNeutron)/MeVtoinvFM << " " << ChemicalPotentialProton/MeVtoinvFM << '\n'; 

		if(BaryonDensity >= 0.999*0.153 && BaryonDensity <= 1.001*0.153){
			cout << BaryonDensity << " " << EnergyDensity << " " << BindingPerNucleon/MeVtoinvFM << " " << '\n';
			print_state(iter, S);				//PRINT THE CURRENT STATUS
			print_state2(iter2, S2);
			cout << "status1 = " << gsl_strerror(status) << '\n';
			cout << "status2 = " << gsl_strerror(status2) << '\n';
			cout << "Asymmetry: " << AlphaAsymmetry << ", Neutron Density: " << NeutronDensity << ", Baryon Density: " << BaryonDensity << 
				", Proton Density: " << ProtonDensity << '\n' << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << '\n';
			cout << '\n' << '\n' << EOmega << " :OMEGA: " << EOmega2 << '\n' <<
				ERho << " :RHO: " << ERho2 << '\n' <<
				ESigma << " :SIGMA: " << ESigma2 << '\n' <<
				ScalarDensity << " :nS: " << ScalarDensity2 << '\n' << 
				Ekin << " :EK: " << Ekin2 << '\n';
		}


		gsl_multiroot_fsolver_free(S2);			//CLEAR MEMORY
		gsl_vector_free(y);				//ON EVERYTHING



		BaryonDensity*= 1.001;
	}
		gsl_multiroot_fsolver_free(S);			//CLEAR MEMORY
		gsl_vector_free(x);				//ON EVERYTHING

	DATER.close();
	return 0;
	
}


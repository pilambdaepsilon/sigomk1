#include "sigointegratorsDM.hh"

int main(){
	double KineticEnergy, KineticEnergyProton, KineticEnergyNeutron, EnergyDensity, Int1, Int2, Int3, IntP1;
	double E0 = 0;	
	double BaryonDensity = 0.001;				//the one true var
	double ConstC, ConstB, GSIGMA, GOMEGA, GRHO;

	double pr1 = 0;
	double pr2 = 0;
	double pr3 = 0;
	double pr4 = 0;
	double pr5 = 0;
	double pr6 = 0;
	double pr7 = 0;
	double pr8 = 0;
	double PARAMSET[6][9] = {{0}};
	string fno;
	cout << '\n' << '\n' << "FILENUMBER: "; cin >> fno;
	ifstream dater;
	dater.open("nchi.dat");
	if (!dater){
		cout << "not open" << '\n';
	}

	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 >> pr4 >> pr5 >> pr6 >> pr7 >> pr8;
	//	cout << pr1 << ' ' << pr2 << ' ' << pr3 << ' ' << pr4 << ' ' << pr5 << ' ' << pr6 << ' ' << pr7 << ' '  << pr8 << '\n';
		if (pr1 == 10){
			PARAMSET[0][0] = pr1;
			PARAMSET[0][1] = pr2;
			PARAMSET[0][2] = pr3;
			PARAMSET[0][3] = pr4;
			PARAMSET[0][4] = pr5;
			PARAMSET[0][5] = pr6;
			PARAMSET[0][6] = pr7;
			PARAMSET[0][7] = pr8;
		}
		else if (pr1 >= 5e9 && pr1 <= 6e9){
			PARAMSET[1][0] = pr1;
			PARAMSET[1][1] = pr2;
			PARAMSET[1][2] = pr3;
			PARAMSET[1][3] = pr4;
			PARAMSET[1][4] = pr5;
			PARAMSET[1][5] = pr6;
			PARAMSET[1][6] = pr7;
			PARAMSET[1][7] = pr8;
		}
		else if (pr1 >= 6e9 && pr1 <= 7e9){
			PARAMSET[2][0] = pr1;
			PARAMSET[2][1] = pr2;
			PARAMSET[2][2] = pr3;
			PARAMSET[2][3] = pr4;
			PARAMSET[2][4] = pr5;
			PARAMSET[2][5] = pr6;
			PARAMSET[2][6] = pr7;
			PARAMSET[2][7] = pr8;
		}
		else if (pr1 >= 7e8 && pr1 <= 8e9){
			PARAMSET[3][0] = pr1;
			PARAMSET[3][1] = pr2;
			PARAMSET[3][2] = pr3;
			PARAMSET[3][3] = pr4;
			PARAMSET[3][4] = pr5;
			PARAMSET[3][5] = pr6;
			PARAMSET[3][6] = pr7;
			PARAMSET[3][7] = pr8;
		}
		else if (pr1 >= 8e9 && pr1 <= 9e9){
			PARAMSET[4][0] = pr1;
			PARAMSET[4][1] = pr2;
			PARAMSET[4][2] = pr3;
			PARAMSET[4][3] = pr4;
			PARAMSET[4][4] = pr5;
			PARAMSET[4][5] = pr6;
			PARAMSET[4][6] = pr7;
			PARAMSET[4][7] = pr8;
		}
		else if (pr1 >= 9e9 && pr1 <= 1e10){
			PARAMSET[5][0] = pr1;
			PARAMSET[5][1] = pr2;
			PARAMSET[5][2] = pr3;
			PARAMSET[5][3] = pr4;
			PARAMSET[5][4] = pr5;
			PARAMSET[5][5] = pr6;
			PARAMSET[5][6] = pr7;
			PARAMSET[5][7] = pr8;
		}
	}


	dater.close();

	ofstream DAT;						//my file
	string sname;
	DAT.open(string("SigOm"+fno+"En.dat").c_str());
	ofstream DATCHECK;
	DATCHECK.open("check.dat");

	double drho = 0.001*rho0;
	double Etest = 0;
	double kFermi, kFermiproton, kFermineutron;
	double ScalarDMDensity = 0;
	double Pressure = 0;
	double kFermi2 = 0;
	double SRtest = 0;
	double ScalarDensity = 0;
	double NeutronDensity = 0;
	double ProtonDensity = 0;
	double BetaRatio = 0;
	double DelMass = 0;
	double EffMass = mass;
	double EffMassPrime = mstar;
	double E01 = mass/MeVtoinvFM;
	int calls = 1000;
	int DMTimeSlice = 5;
	cout << "Beta Ratio (Neutron Density - Proton Density)/Total Density: "; cin >> BetaRatio; cout << '\n';
	cout << "NSage - [0]10 yr. [1]~1e5 yr. [2]~1e8 yr: " ; cin >> DMTimeSlice; 

	double GPI = 0;
	double GCHI = 0;
	double NSage = PARAMSET[DMTimeSlice][0];			//years
	double DMDensity = PARAMSET[DMTimeSlice][2];			//fm^-3
	double SIGCHIN = PARAMSET[DMTimeSlice][3]*1e26;			//fm^2
	double SIGCHI2 = PARAMSET[DMTimeSlice][4]*1e26;			//fm^2
	double DMmass = PARAMSET[DMTimeSlice][5]*1e3*MeVtoinvFM;	//fm^-1

	cout << NSage << " years" << '\n' << '\n';
/*========== now use these to get the variable parameters ==================*/	
	while (BaryonDensity <= 15*rho0){
		kFermi = pow((6*pi*pi/4. *BaryonDensity), (1/3.));
		double kFermiDM = pow((6*pi*pi/4. *DMDensity), (1/3.));
		NeutronDensity = 0.5*(1 + BetaRatio)*BaryonDensity;
		ProtonDensity = 0.5*(1 - BetaRatio)*BaryonDensity;

		kFermiproton = pow((6*pi*pi/4.* ProtonDensity), 1/3.);
		kFermineutron = pow((6*pi*pi/4.* NeutronDensity), 1/3.);

		Int1 = Int_eye1(0, kFermi, calls);
		Int2 = Int_eye2(0, kFermi, calls);
		Int3 = Int_eye3(0, kFermi, calls);

		IntP1 = Int_pee1(0, kFermi, calls);

		ConstC = C(BaryonDensity, Int1, Int2, Int3, kFermi);
		ConstB = B(BaryonDensity, Int1, Int2, Int3, kFermi);
		GSIGMA = Gsigma2(BaryonDensity, ConstB, ConstC, Int1, kFermi);
		GOMEGA = Gomega2(BaryonDensity, kFermi);
		GRHO = Grho2(BaryonDensity, kFermi);
		EnergyDensity = 0;
		for(int i = 0; i < 100; i++){
			EffMassPrime = EffMass;
			ScalarDensity = Int_rhoS(0, kFermi, EffMassPrime, calls);
			EffMass = mass - GSIGMA * ScalarDensity;
		}

		Int1 = Int_eye1(0, kFermi, calls);
		Int2 = Int_eye2(0, kFermi, calls);
		Int3 = Int_eye3(0, kFermi, calls);
		IntP1 = Int_pee1(0, kFermi, calls);

		ConstC = C(BaryonDensity, Int1, Int2, Int3, kFermi);
		ConstB = B(BaryonDensity, Int1, Int2, Int3, kFermi);
		GSIGMA = Gsigma2(BaryonDensity, ConstB, ConstC, Int1, kFermi);
		GOMEGA = Gomega2(BaryonDensity, kFermi);
		GRHO = Grho2(BaryonDensity, kFermi);
		
		GCHI = Gchi(SIGCHI2, DMmass, MPI);
		GPI = Gpi(SIGCHIN, mass, SIGCHI2, MPI, GCHI);
		DelMass = GSIGMA*ScalarDensity;
		KineticEnergy = Int_Ekinetic(0, kFermi, EffMass, calls);
		KineticEnergyProton = Int_Ekinetic(0, kFermiproton, EffMass, calls);
		KineticEnergyNeutron = Int_Ekinetic(0, kFermineutron, EffMass, calls);
		double KineticEnergyDM = Int_Ekinetic(0, kFermiDM, DMmass, calls);

		double EnergyDensitynot = 
		 KineticEnergy
		+ 0.5* GOMEGA * pow(BaryonDensity, 2)		
		+ (1./8.)* GRHO * pow((ProtonDensity - NeutronDensity), 2) 
		+ 0.5* GSIGMA * pow(ScalarDensity, 2)
		+ ( (1./3.)*ConstB*mass*pow(DelMass, 3) + (0.25) * ConstC*pow(DelMass, 4) );

		if (BaryonDensity <= 0.001 && BaryonDensity > 0.0){
			E01 = EnergyDensitynot/BaryonDensity/MeVtoinvFM;
			//E01 = mass/MeVtoinvFM;
			cout << "CHECK" << '\n' << '\n';
		}

		double EnergyDensity = 
		 KineticEnergy
		+ 0.5* GOMEGA * pow(BaryonDensity, 2)		
		+ (1./8.)* GRHO * pow((ProtonDensity - NeutronDensity), 2) 
		+ 0.5* GSIGMA * pow(ScalarDensity, 2)
		+ ( (1./3.)*ConstB*mass*pow(DelMass, 3) + (0.25) * ConstC*pow(DelMass, 4) );

		EnergyDensity += 0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
		+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
		+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity;
		+ KineticEnergyDM; 

		double Pressure =
		(1./3.)* IntP1 +
		+ 0.5* GOMEGA * pow(BaryonDensity, 2)		
		+ (1./8.)* GRHO * pow((ProtonDensity - NeutronDensity), 2) 
		- 0.5* GSIGMA * pow(ScalarDensity, 2)
		- ( (1./3.)*ConstB*mass*pow(DelMass, 3) - (0.25) * ConstC*pow(DelMass, 4) );


		double BindingPerNucleon = EnergyDensity/(BaryonDensity*MeVtoinvFM) - E01;


		/* Assuming an attractive potential between DM particles */
		Pressure -= (0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
		+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
		+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity);

		BindingPerNucleon *=6;

		DAT <<							//0 
		BaryonDensity << " " <<					//1
		BindingPerNucleon << " " <<			//2
		(-Pressure/MeVtoinvFM) << '\n'; //" " <<

/*		EnergyDensity/MeVtoinvFM << " " <<
		GSIGMA << " " <<					//3
		GOMEGA << " " <<					//4
		GRHO << " " <<						//5
		GSIGMA*ScalarDensity << " " << 				//6
		GOMEGA*BaryonDensity << " " << 				//7
		ScalarDensity << " " << 				//8
		EffMass/mass << " " <<					//9
		mass << " " <<						//10
		kFermi << '\n';						//11
		
*/
		if (BaryonDensity < 1.00*rho0 && BaryonDensity > 0.999*rho0){
			cout << 
			"r/ro: " << BaryonDensity/rho0<< '\n' << 
			"Scalar Density: " << ScalarDensity << '\n' << 
			"Neutron/ Proton Densities: " << NeutronDensity << "/ " << ProtonDensity << '\n' << 
			"EnergyDensity: " << EnergyDensity  << '\n' << 
			"E_the: " << (BperA + mass)*rho0 << '\n' << '\n' <<
			"BindingPerNucleon: " << BindingPerNucleon << '\n' << '\n' <<
			"K: " << Kompress/MeVtoinvFM << '\n' << 
			"sig: " << GSIGMA << '\n' <<
			"om: " << GOMEGA << '\n' <<
			"rho: " << GRHO << '\n' <<
			"b: " << ConstB << '\n' << 
			"c: " << ConstC << '\n' << 
			"kFermi(P): " << kFermiproton << '\n' <<  
			"kFermi(N): " << kFermineutron << '\n' <<  
			"kFermi: " << kFermi << '\n' <<  
			"Ek: " << KineticEnergy << '\n' <<
			"Eo: " << E01 << '\n' <<
			"m: " << mass/MeVtoinvFM << "[MeV], " << mass << "[fm^-1]" << "        ______" << '\n' <<
			"m*: " << EffMass/MeVtoinvFM << "[MeV], " << EffMass << "[fm^-1]" << " _____|---> " << mstar/mass << '\n' <<
			"Sigchin: " << PARAMSET[DMTimeSlice][3] << " [cm^2]" << '\n' << 
			"Sigchi2: " << PARAMSET[DMTimeSlice][4] << " [cm^2]" << '\n' <<
			"DMmass: " << DMmass/1e3/MeVtoinvFM << " [GeV]" << '\n' << '\n';
			
			
			double i1 = 0;
			i1 = eye3(kFermi, EffMass);
				
			cout << "int: " << Int3 << "    an: " << i1 << '\n' << '\n';
		}
		
		BaryonDensity += drho;
	}

	DAT.close();
	return 0;
} 

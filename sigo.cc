#include "sigointegrators.hh"

int main(){
	double KineticEnergy, KineticEnergyProton, KineticEnergyNeutron, EnergyDensity, Int1, Int2, Int3;
	double E0 = 0;	
	double BaryonDensity = 0.001;				//the one true var
	double ConstC, ConstB, GSIGMA, GOMEGA, GRHO;

	ofstream DAT;						//my file
	DAT.open("SigOm.dat");
	ofstream DATCHECK;
	DATCHECK.open("check.dat");

	double drho = 0.001*rho0;
	double Etest = 0;
	double kFermi, kFermiproton, kFermineutron;
	double ScalarDMDensity = 0;
	double DMDensity = 0;
	double GCHI, GPI, DMmass, PImass;
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
	cout << "Beta Ratio (Neutron Density - Proton Density)/Total Density: "; cin >> BetaRatio; cout << '\n';

/*========== now use these to get the variable parameters ==================*/	
	while (BaryonDensity <= 10*rho0){
		kFermi = pow((6*pi*pi/4. *BaryonDensity), (1/3.));
		NeutronDensity = 0.5*(1 + BetaRatio)*BaryonDensity;
		ProtonDensity = 0.5*(1 - BetaRatio)*BaryonDensity;

		kFermiproton = pow((6*pi*pi/4.* ProtonDensity), 1/3.);
		kFermineutron = pow((6*pi*pi/4.* NeutronDensity), 1/3.);

		Int1 = Int_eye1(0, kFermi, calls);
		Int2 = Int_eye2(0, kFermi, calls);
		Int3 = Int_eye3(0, kFermi, calls);

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
			
			// CHECK FOR CONVERGENCE OF EFFECTIVE MASS AND SCALAR DENSITY
			if(DATCHECK.is_open()){
				DATCHECK << i << " " << EffMassPrime/mass << " " << ScalarDensity << '\n';	
			}
		}
		DATCHECK.close();
		
		DelMass = GSIGMA*ScalarDensity;
		KineticEnergy = Int_Ekinetic(0, kFermi, EffMass, calls);
		KineticEnergyProton = Int_Ekinetic(0, kFermiproton, EffMass, calls);
		KineticEnergyNeutron = Int_Ekinetic(0, kFermineutron, EffMass, calls);

		EnergyDensity = 
		- KineticEnergy
		- 0.5* GOMEGA * pow(BaryonDensity, 2)		
		- (1./8.)* GRHO * pow((ProtonDensity - NeutronDensity), 2) 
		+ 0.5* GSIGMA * pow(ScalarDensity, 2)
		+ ( (1./3.)*ConstB*mass*pow(DelMass, 3) + (0.25) * ConstC*pow(DelMass, 4) );


		if (BaryonDensity <= 0.001 && BaryonDensity > 0.0){
			E01 = EnergyDensity/BaryonDensity/MeVtoinvFM;
			//E01 = mass/MeVtoinvFM;
			cout << "CHECK" << '\n' << '\n';
		}
	
		double BindingPerNucleon = EnergyDensity/(BaryonDensity*MeVtoinvFM) - E01;

		DAT <<							//0 
		BaryonDensity << " " <<					//1
		BindingPerNucleon << " " <<				//2
		GSIGMA << " " <<					//3
		GOMEGA << " " <<					//4
		GRHO << " " <<						//5
		GSIGMA*ScalarDensity << " " << 				//6
		GOMEGA*BaryonDensity << " " << 				//7
		ScalarDensity << " " << 				//8
		EffMass/mass << " " <<					//9
		kFermi << '\n';						//10
		
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
			"m*: " << EffMass/MeVtoinvFM << "[MeV], " << EffMass << "[fm^-1]" << " _____|---> " << mstar/mass << '\n' << '\n';
			
			double i1 = 0;
			i1 = eye3(kFermi, EffMass);
				
			cout << "int: " << Int3 << "    an: " << i1 << '\n' << '\n';
		}
		
		BaryonDensity += drho;
	}

	DAT.close();
	return 0;
} 

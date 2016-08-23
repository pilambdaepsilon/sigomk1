#include "sigointegrators.hh"

int main(){
	double KineticEnergy, EnergyDensity, Int1, Int2, Int3;			
	double E0 = 0;	
	double BaryonDensity = 0.0001;				//the one true variable
	double ConstC, ConstB, GSIGMA, GOMEGA, GRHO;

	ofstream DAT;						//write data into this file
	DAT.open("SigOm.dat");
	ofstream DATCHECK;
	DATCHECK.open("check.dat");

	double drho = 0.001*rho0;
	double kFermi = 0;
	double ScalarDMDensity, DMDensity, GCHI, GPI, DMmass, PImass;
	double ScalarDensity, NeutronDensity, ProtonDensity, BetaRatio, DelMass;
	double EffMass = 0;
	double EffMassPrime = mass;
	double E01 = mass/MeVtoinvFM;
	int calls = 1000;

	cout << "Beta Ratio (Neutron Density - Proton Density)/Total Density: "; cin >> BetaRatio; cout << '\n';

/*========== now use these to get the variable parameters ==================*/	
	while (BaryonDensity <= 10*rho0){
		kFermi = pow((6*pi*pi/4 *BaryonDensity), (1./3));
		Int1 = Int_eye1(0, kFermi, calls);
		Int2 = Int_eye2(0, kFermi, calls);
		Int3 = Int_eye3(0, kFermi, calls);

		ConstC = C(BaryonDensity, Int1, Int2, Int3, kFermi);
		ConstB = B(BaryonDensity, Int1, Int2, Int3, kFermi);
		GSIGMA = Gsigma2(BaryonDensity, ConstB, ConstC, Int1, kFermi);
		GOMEGA = Gomega2(BaryonDensity, kFermi);
		EnergyDensity = 0;

		for(int i = 0; i < 25; i++){
			EffMass = EffMassPrime;
			ScalarDensity = Int_rhoS(0, kFermi, EffMass, calls);
			EffMassPrime = mass - GSIGMA * ScalarDensity;
			
			// CHECK FOR CONVERGENCE OF EFFECTIVE MASS AND SCALAR DENSITY
			if(DATCHECK.is_open()){
				DATCHECK << i << " " << EffMassPrime/mass << " " << ScalarDensity << '\n';	
			}
		}
		DATCHECK.close();
		
		GRHO = Grho2(BaryonDensity, kFermi);
		DelMass = mass - EffMassPrime;

		NeutronDensity = 0.5*(1+BetaRatio)*BaryonDensity;
		ProtonDensity = 0.5*(1 - BetaRatio)*BaryonDensity;

		double kFermiproton = pow((6.*pi*pi/4.* ProtonDensity), 1./3);
		double kFermineutron = pow((6.*pi*pi/4.* NeutronDensity), 1./3);
		double KineticEnergyProton = Int_Ekinetic(0, kFermiproton, calls);
		double KineticEnergyNeutron = Int_Ekinetic(0, kFermineutron, calls);

		EnergyDensity = 
		- KineticEnergyProton/(2.*BaryonDensity)
		- KineticEnergyNeutron/(2.*BaryonDensity)
		- 0.5* GOMEGA * BaryonDensity		
		- 0.5* GRHO * pow((ProtonDensity - NeutronDensity), 2.)/BaryonDensity 
		- 0.5*(GSIGMA)*pow(ScalarDensity, 2.)/BaryonDensity 
		- (1./3.)*ConstB*mass*pow(DelMass, 3.)/BaryonDensity - (1./4.) * ConstC*pow(DelMass, 4.)/BaryonDensity;


		if (BaryonDensity <= 0.001 && BaryonDensity > 0.0){
			E01 = EnergyDensity/MeVtoinvFM;
		}
	
		DAT <<							//0 
		BaryonDensity<< " " << 					//1
		(EnergyDensity)/MeVtoinvFM - E01 << " " << 		//2
		GSIGMA << " " <<					//3
		GOMEGA << " " <<					//4
		GRHO << " " <<						//5
		GSIGMA*ScalarDensity<< " " << 				//6
		-GOMEGA*BaryonDensity<< " " << 				//7
		ScalarDensity << " " << 				//8
		EffMassPrime/mass << " " <<				//9
		kFermi << '\n';						//10
		
		if (BaryonDensity < 1.00*rho0 && BaryonDensity > 0.999*rho0){
			cout << 
			"r/ro: " << BaryonDensity/rho0<< '\n' << 
			"E: " <<(EnergyDensity)/MeVtoinvFM - E01  << '\n' << 
			"K: " << Kompress/MeVtoinvFM << '\n' << 
			"sig: " << GSIGMA << '\n' <<
			"om: " << GOMEGA << '\n' <<
			"rho: " << GRHO << '\n' <<
			"b: " << ConstB << '\n' << 
			"c: " << ConstC << '\n' << 
			"kFermi(P): " << kFermiproton << '\n' <<  
			"kFermi(N): " << kFermineutron << '\n' <<  
			"kFermi: " << kFermi << '\n' <<  
			"Eo: " << E01 << '\n' <<
			"m: " << mass/MeVtoinvFM << "    ______" << '\n' <<
			"m*: " << mstar/MeVtoinvFM << " _____|---> " << mstar/mass << '\n' << '\n';
			
			double i1 = 0;
			i1 = eye3(kFermi, mstar);
				
			cout << "int: " << Int3 << "    an: " << i1 << '\n' << '\n';
		}
		
		BaryonDensity += drho;
	}

	DAT.close();
	return 0;
} 

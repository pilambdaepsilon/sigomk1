#include "sigointegrators.hh"

int main(){
	double ScalarDensity, KineticEnergy, EnergyDensity, Int1, Int2, Int3;			
	double E0 = 0;	
	double BaryonDensity = 0;				//the one true var
	double ConstC, ConstB, GSIGMA, GOMEGA, GRHO;
	ofstream DAT;						//my file
	DAT.open("SigOm.dat");

	double drho = 0.01*rho0;
	double Etest = 0;
	double kFermi = 0;
	double kFermi2 = 0;
/*========== now use these to get the variable parameters ==================*/	
	while (BaryonDensity <= 10*rho0){
		kFermi = pow((6*pi*pi/4 *BaryonDensity), 0.3333333);
		ScalarDensity = Int_rhoS(0, kFermi, 1000);
		KineticEnergy = Int_Ekinetic(0, kFermi, 1000);
		Int1 = Int_eye1(0, kFermi, 1000);
		Int2 = Int_eye2(0, kFermi, 1000);
		Int3 = KineticEnergy;

		ConstC = B(BaryonDensity, Int1, Int2, Int3, kFermi);
		ConstB = C(BaryonDensity, Int1, Int2, Int3, kFermi);
		GSIGMA = Gsigma2(BaryonDensity, ConstB, ConstC, Int1, kFermi);
		GOMEGA = Gomega2(BaryonDensity, kFermi);
		GRHO = Grho2(BaryonDensity, kFermi);
		EnergyDensity = 0;
		
//		EnergyDensity += 0.5 * GSIGMA * pow(ScalarDensity,2);
//		EnergyDensity += 0.5 * GOMEGA * pow(BaryonDensity,2);
//		EnergyDensity += 0.5 * GRHO * pow(BaryonDensity,2);
		EnergyDensity += KineticEnergy;
//		EnergyDensity += 0.3333* ConstB*mass* pow(GSIGMA, 3)*pow(ScalarDensity, 3);
//		EnergyDensity += 0.25* ConstC * pow(GSIGMA, 4)*pow(ScalarDensity,4);

		Etest = 3*kFermi*kFermi/(10*mass);
//		Etest +=0.5*(GOMEGA - GSIGMA)*BaryonDensity;
//		Etest += GSIGMA*BaryonDensity/mass *(3*pow(kFermi, 3)/(10*mass));

		if (BaryonDensity/rho0 < 1.01 && BaryonDensity/rho0 > 0.99){
			cout << 
			"r/ro: " << BaryonDensity/rho0<< '\n' << 
			"E: " << (EnergyDensity/BaryonDensity/MeVtoinvFM)  << '\n' << 
			"Etest: " << Etest/MeVtoinvFM << '\n' <<
			"sig: " << GSIGMA << '\n' <<
			"om: " << GOMEGA << '\n' <<
			"rho: " << GRHO << '\n' << 
			"b: " << B(BaryonDensity, Int1, Int2, Int3, kFermi) << '\n' << 
			"c: " << C(BaryonDensity, Int1, Int2, Int3, kFermi) << '\n' << 
			"b2: " << B2(BaryonDensity, Int1, Int2, Int3, kFermi) << '\n' << 
			"c2: " << C2(BaryonDensity, Int1, Int2, Int3, kFermi) << '\n' << '\n';
		}

		DAT << 
		BaryonDensity/rho0<< " " << 
		(EnergyDensity)  << " " << 
		Etest << " " <<
		GSIGMA << " " <<
		GOMEGA << " " <<
		GRHO << " " <<
		ConstB << " " <<
		ConstC << " " <<
		ScalarDensity << " " << 
		kFermi << '\n'; 
		
		BaryonDensity += drho;
	}

	DAT.close();
	return 0;
} 



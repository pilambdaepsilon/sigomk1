#include "sigointegrators.hh"

int main(){
	double KineticEnergy, EnergyDensity, Int1, Int2, Int3;			
	double E0 = 0;	
	double BaryonDensity = 0*rho0;				//the one true var
	double ConstC, ConstB, GSIGMA, GOMEGA, GRHO;
	ofstream DAT;						//my file
	DAT.open("SigOm.dat");

	double drho = 0.01*rho0;
	double Etest = 0;
	double kFermi = 0;
	double kFermi2 = 0;
	double SRtest = 0;
	double ScalarDensity = 0;
	double DelMass = 0;
	double E01 = mass/MeVtoinvFM;
	double E02 = mass/MeVtoinvFM;
	double supress = 0;
	cout << "supress rho: "; cin >> supress; cout << '\n';
/*========== now use these to get the variable parameters ==================*/	
	while (BaryonDensity <= 10*rho0){
		kFermi = pow((6*pi*pi/4 *BaryonDensity), (1./3));
		KineticEnergy = Int_Ekinetic(0, kFermi, 1000);
		Int1 = Int_eye1(0, kFermi, 1000);
		Int2 = Int_eye2(0, kFermi, 1000);
		Int3 = KineticEnergy;

		ConstC = C(BaryonDensity, Int1, Int2, Int3, kFermi);
		ConstB = B(BaryonDensity, Int1, Int2, Int3, kFermi);
		GSIGMA = Gsigma2(BaryonDensity, ConstB, ConstC, Int1, kFermi);
		GOMEGA = Gomega2(BaryonDensity, kFermi);
		EnergyDensity = 0;
		Etest = 0;
 		SRtest = BaryonDensity *(1 - 3*kFermi*kFermi/(10*mass*mass));
		ScalarDensity = Int_rhoS(0, kFermi, 1000);
		DelMass = mass - mstar;
		GRHO = Grho2(BaryonDensity, kFermi);
	
		EnergyDensity = 0.5*(GOMEGA)* BaryonDensity - 0.5*(GSIGMA)*pow(ScalarDensity, 2)/BaryonDensity + KineticEnergy/BaryonDensity - (1./3)*ConstB*mass*pow((mass-mstar), 3)/BaryonDensity - (1./4) * ConstC*pow(DelMass, 4)/BaryonDensity
		+ 0.5*GRHO*BaryonDensity*supress; 
		Etest = 0.5*(GOMEGA)* BaryonDensity - 0.5*(GSIGMA)*pow(SRtest, 2)/BaryonDensity + KineticEnergy/BaryonDensity - 0.3333*ConstB*mass*pow((mass-mstar), 3)/BaryonDensity - 0.25 * ConstC*pow(DelMass, 4)/BaryonDensity
	+ 0.5*GRHO*BaryonDensity*supress;

		if (BaryonDensity < 0.02 && BaryonDensity > 0.0){
			E01 = EnergyDensity/MeVtoinvFM;
			E02 = Etest/MeVtoinvFM;
		}
	
		DAT << 
		BaryonDensity/rho0<< " " << 
		(EnergyDensity)/MeVtoinvFM - E01 << " " << 
		(Etest)/MeVtoinvFM - E02<< " " <<
		-GSIGMA << " " <<
		GOMEGA << " " <<
		GRHO << " " <<
//		ConstB << " " <<
//		ConstC << " " <<
		GSIGMA*ScalarDensity*ScalarDensity<< " " << 
		Int3 << " " << '\n';
		if (BaryonDensity < 1.01*rho0 && BaryonDensity > 0.99*rho0){
			cout << 
			"r/ro: " << BaryonDensity/rho0<< '\n' << 
			"E: " <<(EnergyDensity)/MeVtoinvFM - E01  << '\n' << 
			"Etest: " << Etest/MeVtoinvFM - E02 << '\n' <<
			"K: " << Kompress/MeVtoinvFM << '\n' << 
			"sig: " << GSIGMA << '\n' <<
			"om: " << GOMEGA << '\n' <<
			"rho: " << GRHO << '\n' <<
			"b: " << ConstB << '\n' << 
			"c: " << ConstC << '\n' << 
			"kFermi: " << kFermi << '\n' << '\n';
			
			double i1 = 0;
			i1 = eye3(kFermi, mstar);
				
			cout << "int: " << Int3 << "    an: " << i1 << '\n' << '\n';
		}
		
		BaryonDensity += drho;
	}

	DAT.close();
	return 0;
} 



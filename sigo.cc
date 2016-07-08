#include "sigointegrators.hh"

int main(){
	double ScalarDensity = 0;
	double BaryonDensity = 0;
	double KineticEnergy = 0;
	double EnergyDensity = 0;
	double E0 = 0;
	ofstream DAT;
	DAT.open("SigOm.dat");
	ScalarDensity = Int_rhoS(0, kFermi, 1000);
	KineticEnergy = Int_Ekinetic(0, kFermi, 1000);

/*****
       /$$           /$$                                     /$$                    
      | $$          | $$                                    |__/                    
  /$$$$$$$  /$$$$$$ | $$$$$$$  /$$   /$$  /$$$$$$   /$$$$$$  /$$ /$$$$$$$   /$$$$$$ 
 /$$__  $$ /$$__  $$| $$__  $$| $$  | $$ /$$__  $$ /$$__  $$| $$| $$__  $$ /$$__  $$
| $$  | $$| $$$$$$$$| $$  \ $$| $$  | $$| $$  \ $$| $$  \ $$| $$| $$  \ $$| $$  \ $$
| $$  | $$| $$_____/| $$  | $$| $$  | $$| $$  | $$| $$  | $$| $$| $$  | $$| $$  | $$
|  $$$$$$$|  $$$$$$$| $$$$$$$/|  $$$$$$/|  $$$$$$$|  $$$$$$$| $$| $$  | $$|  $$$$$$$
 \_______/ \_______/|_______/  \______/  \____  $$ \____  $$|__/|__/  |__/ \____  $$
                                         /$$  \ $$ /$$  \ $$               /$$  \ $$
                                        |  $$$$$$/|  $$$$$$/              |  $$$$$$/
                                         \______/  \______/                \______/ 	****/


	double Ssig = 0;
	double ctest = 0;
	double btest = 0;
	cout << "Field Strength: "; cin >> Ssig;
	cout << "c: "; cin >> ctest;
	cout << "b: "; cin >> btest;
	
	double ratio = 1/Ssig + 3*pow((mass - mstar),2)*ctest + 2*pow((mass-mstar),2)*btest;
	double ratio2 = del1(rho0)/alpha1(rho0);
	double eta = 6*pow(kFermi, 3)/(pi*pi) * mstar*mstar/(mstar*mstar + kFermi2);

	double Alf = (-I1 - 1/Ssig - 3*pow((mass - mstar),2)*ctest - 2*pow((mass-mstar),2)*btest);
 
	
	cout << '\n' << "exp ratio: " << ratio << "   ac ratio: " << ratio2 << '\n';
	cout << "alpha1: " << alpha1(rho0) << '\n';
	cout << "delta1: " << del1(rho0) << '\n';
	cout << "Alf: " << Alf << '\n';
	cout << "I1: " << I1 << '\n';


/*======================================================================================================*/

	while (BaryonDensity <= 10*rho0){
		EnergyDensity = 0;
			
		EnergyDensity = 0.5 * Gsigma2(BaryonDensity) * pow(ScalarDensity,2);
		EnergyDensity += 0.5 * Gomega2(BaryonDensity) * pow(BaryonDensity,2);
		EnergyDensity += 0.5 * Grho2(BaryonDensity) * pow(BaryonDensity,2);
		EnergyDensity += KineticEnergy;
		EnergyDensity += 0.3*B(BaryonDensity)*mass* pow(Gsigma2(BaryonDensity), 3)*pow(ScalarDensity, 3);
		EnergyDensity += 0.25*C(BaryonDensity) * pow(Gsigma2(BaryonDensity), 4)*pow(ScalarDensity,4);

//		if(BaryonDensity == 0.){
//			EnergyDensity = 208.18*MeVtoinvFM;
//		}
		if(BaryonDensity == rho0){
			cout << '\n' << "rho/rho0: " << BaryonDensity/rho0 << '\n' << "EnergyDensity: " << EnergyDensity/MeVtoinvFM - 208.18<< '\n' << "Sigma: " << Gsigma2(BaryonDensity)<< '\n' << "Omega: " << Gomega2(BaryonDensity)<< '\n' << "Rho: " << Grho2(BaryonDensity)<< '\n' << "b: " << B(BaryonDensity)<< '\n' << "c: " << C(BaryonDensity) << '\n' << '\n';
		}
		//208.18
		DAT << BaryonDensity/rho0 << " " << EnergyDensity/BaryonDensity << '\n';
		
		BaryonDensity += rho0*0.25;
	}

	DAT.close();
	return 0;
} 



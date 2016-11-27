#include "dmnofunc.hh"

int main(){
/*================== Pull in Cross-Sections and masses ====================*/
	double DMmass = 0;
	double powsigchin =0;
	double powsigchi2 =0;
	cout << '\n';
	cout << "DM-Nucleon cross-section[order, cm^2]: "; cin >> powsigchin;
	cout << "DM-DM cross-section[order, cm^2]: "; cin >> powsigchi2;
	double SIGCHIN = pow(10., powsigchin);
	double SIGCHI2 = pow(10., powsigchi2);
	cout << "DM mass[GeV]: "; cin >> DMmass;
	double DMmass0 = DMmass;

	if (SIGCHI2/DMmass >= 2e-24){
		while(SIGCHI2/DMmass >=  2e-24){
			cout << '\n' << "Combination Breaks Bullet Cluster... Pick new ones" << '\n';
			cout << "DM-DM cross-section[order, cm^2]: "; cin >> powsigchi2;
			SIGCHI2 = pow(10., powsigchi2);
			cout << "DM mass[GeV]: "; cin >> DMmass;
		}
	}
	double ReducedMass = DMmass * Nmass/(DMmass + Nmass);
	double CaptureRate = 9.19e22 * (SIGCHIN/1e-55) /DMmass * rhoDM/pow(convinvGeVtocm, 3.);				//Capture rate using these parameters
	double SelfCapture = 1.06e-3 * (SIGCHI2/1e-24) /DMmass * rhoDM/pow(convinvGeVtocm, 3.);				//Self-capture rate using these parameters

/*======================= Core Density ========================================*/
	double nbee = 3*MNS/(Nmass*4*pi*pow(RNS,3));
	cout << "nB(pF): " << nbee << " ... nB:" << CoreDens << '\n';
	double TIME = 1;
	ofstream myfile;
	ofstream myfile2;
	ofstream myfile3;
	myfile.open("nchi.dat");
	myfile2.open("Eboundcheck.dat");
	myfile3.open("geotime.dat");
	
	double DMNO = 0;
	double DMNOplus = 0;
	
/*=========================== Calculate the Chandrasekhar Limits ==========================================*/
	double chandpref = pow( (9./(32.*pi*pi)), 1./3.)*Planckbar*clight*5./3. /GNewton;
	double Nchichand = pow( (chandpref/(DMmass*DMmass)), 3./2.);
	double Nchichandf = 1.8e51*pow((100/DMmass), 3.);
	double Nchichandb = 1.5e34*pow((100/DMmass), 2.);

/*======================== Get the TIMESCALES ====================================================*/
	double TCONTAIN = tcon(DMmass, SIGCHIN);
	if (TCONTAIN > 3.5){
		TCONTAIN = 3.5;
	}
	double TGEO = tgeo(DMmass, ReducedMass, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
	double TTHERM = ttherm(DMmass, Nmass, ReducedMass, SIGCHIN);
	double times[4] = {TGEO, TTHERM, TCONTAIN, TMAX};
	double t1, t2, t3 ,t4;
	int counter = 0;

	vector<double> timevec (times, times+4);
	double RTH = radius(DMmass, ReducedMass, TTHERM, SIGCHIN)/convcmtoinvGeV;
	double RTHEST = sqrt(9.* TNS/(4.*pi*GNewton0*CoreDens*Nmass*DMmass))/convcmtoinvGeV;
	cout << "containment time: " << TCONTAIN << '\n'
		<< "thermalization time: " << TTHERM << '\n'
		<< "goemetric time: " << TGEO << '\n'	
		<< "thermal. radius 1: " << RTH << '\n'	
		<< "thermal. radius 2: " << RTHEST << '\n';	
	
/*================= Sort the time scales for proper evolution =====================================*/
	stable_sort(timevec.begin(), timevec.end());
	t1 = timevec[0];
	t2 = timevec[1];
	t3 = timevec[2];
	t4 = timevec[3];
	cout << "t1: " << t1 << '\n'
		<< "t2: " << t2 << '\n'
		<< "t3: " << t3 << '\n'	
		<< "t4: " << t4 << '\n';	
	double tprev = 0;
	int count = 0;

/*================== Calculate the number of particles by regime ================================*/
	while(TIME <= TMAX){
		if (TIME < t1){
			DMNOplus = 0;
		}
 
		else if (TIME >= t1 && TIME < t2){
			if(t1 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 0){
				cout << "T1: Geometric" << '\n';
				count = 1;}
			}
			else if(t1 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 0){
				cout << "T1: Thermalization" << '\n';
				count = 1;}
			}
			else if(t1 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t1 = 10^10 years." << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 0){
				cout << "T1: Regular Evolution" << '\n';
				count = 1;}
			}
		}

		else if (TIME >= t2 && TIME < t3){
			if(t2 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 1){
				cout << "T2: Geometric" << '\n';
				count = 2;}
			}
			else if(t2 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 1){
				cout << "T2: Thermalization" << '\n';
				count = 2;}
			}
			else if(t2 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t2 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 1){
				cout << "T2: Regular Evolution" << '\n';
				count = 2;}
			}
		}

		else if (TIME >= t3 && TIME <t4){
			if(t3 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 2){
				cout << "T3: Geometric" << '\n';
				count = 3;}
			}
			else if(t3 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture);
				if(count == 2){
				cout << "T3: Thermalization" << '\n';
				count = 3;}
			}
			else if(t3 == TMAX){
				DMNOplus = 0;
				cout << "DONE, t3 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if(count == 2){
				cout << "T3: Regular Evolution" << '\n';
				count = 3;}
			}
		}

		else if (TIME >= t4){
			DMNOplus = 0;
			if(counter == 0){
				cout << "DONE, t4 = 10^10 years" << '\n';
				counter++;
			}
		}

	
		DMNO += DMNOplus;
		double Rchi = radius(DMmass, ReducedMass, TIME, SIGCHIN)/convcmtoinvGeV*1e13;		//DM sphere radius in fm
		double Volumechi = 4./3. *pi * pow(Rchi, 3.);
		double DMDENS = DMNO/(Volumechi);			//DM number density in inverse fm^3
		myfile << TIME << " " << DMNO << " " << DMDENS << " " << SIGCHIN << " " << SIGCHI2 << " " << DMmass << " " << Nchichandf << " " << Nchichandb << '\n';
		tprev = TIME;
		TIME *= 1.01;
	}

/*================================================= OTHER TESTS =====================================================*/
	while(DMmass <= 1e4){
		double Einf = 2*Nmass*DMmass/pow((DMmass + Nmass), 2) * ( 1 - sqrt(1 - 2*GNewton*MNS/(RNS*clight*clight)));
		double Emin = 1./3.*v2/(clight*clight);
		myfile2 << DMmass << " " << Emin << " " << Einf << '\n';
		DMmass *= 10;
	}

	TIME = 3.5; 
	ReducedMass = DMmass0 * Nmass/(DMmass0 + Nmass);
	while(TIME <= 1e10){
		double rad = radius(DMmass0, ReducedMass, TIME, SIGCHIN)/RNS;
		double rfid = sqrt(Nchi0(CaptureRate, SelfCapture, TIME)* SIGCHI2*pow(convcmtoinvGeV,2.)/(pi * RNS *RNS));

		myfile3 << TIME-3.5 << " " << rad << " " << rfid << '\n';	
//		cout << DMmass0 << ' ' << SIGCHIN << '\n';
		TIME *= 1.001;

	}
	myfile3.close();
	myfile2.close();
	myfile.close();
	return 0;
}

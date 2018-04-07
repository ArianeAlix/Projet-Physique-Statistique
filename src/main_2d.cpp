#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
using namespace std;


//Param�tres de la simulation
const int N = 10;//->nombre de particules = N^2
const double T = 120;//temp�rature du syst�me
const double k = 1.3806485279*pow(10, -23);//constante de boltzmann

//Param�tres des particules
const double sigR = 3.405*pow(10, -10);//sigma r�elle
const double sig = 1;//sigma r�duit
const double a = sig;//taille d'une cellules en sigma
const double epsR = k * T; //vrai epsilon
const double eps = 1;//eps r�duit
const double mR = 6.63*pow(10, -26);//masse r�elle d'un atome d'argon
const double m = 1;//masse r�duite
const double L = N * sig* a; //->c�t� d'un carr�

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.5*sig;

// Constantes li�es au temps
const double deltaT = pow(10, -2);//pas de temps de la simulation
const double tau = deltaT * 10;//pas d'�chantillonnage des positions
const double tmax = 100; //temps de la simulation

void InitRandom()
{
	srand((unsigned int)time(0));
}

double Random(double a, double b)
{
	double x = double(rand()) / RAND_MAX;
	return a + (b - a)*x;
}


double *dist2D(double* qi, double* qj) {
	double dist[2];
	//Doit prendre un compte la pr�iodicit�
	dist[0] = qj[0] - qi[0] - L * round(((qj[0] - qi[0]) / L));
	dist[1] = qj[1] - qi[1] - L * round(((qj[1] - qi[1]) / L));
	return dist;
}


void init2D(double q0[N][2] , double q1[N][2], double v[N][2]) {
	//On choisit un segment de longueur N*a*sigma
	//On place les particules de mani�re r�guli�re sur le segment
	//Et de fa�on � ce que la p�riodicit� implique une r�gularit� 
	//des distances entre toutes les particules
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
	double en = 0.0;//energy of the whole system
	double sumv[2] = { 0,0 };
	double sumvquad = 0;
	for (int i = 0; i < N; i++) {
		q1[i][0] = i * L/N + L/(2.0*N);
		q1[i][1] = i * L / N + L / (2.0*N);
		v[i][0] = distribution(generator);
		v[i][0] = distribution(generator);
		sumv[0] = +v[i][0];
		sumv[1] = +v[i][1];
		sumvquad += pow(v[i][0], 2) + pow(v[i][1], 2);
	}
	double moyv[2] = { sumv[0] / N,sumv[0] / N };
	double moyvquad = sumvquad / N;

	for (int i = 0; i < N; i++) {
		v[i][0] += -moyv[0]; 
		v[i][1] += -moyv[1]; //On fixe la moyenne � 0

		q0[i][0] = q1[i][0] - v[i][0] * deltaT;
		q0[i][1] = q1[i][1] - v[i][1] * deltaT;
	}
}

double pot(double rij) {
	if (rij <= rm) {
		double p = 4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6));
		return p;
	}
	else if (rm < rij && rij  <= rcut) {
		double p = 4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6));
		double x = (rij - rm) / (rcut - rm);
		p *= (2 * pow(x, 3) - 3 * pow(x, 2) + 1);
		return p;
	}
	return 0.0;

}

double pot_prime(double rij) {
	double h = pow(10, -15);
	if (rij <= rcut) {
		double pp = (pot(rij + h) - pot(rij-h)) / (2*h);
		return pp;
	}
	return 0.0;
}


double f1D(double rij) {
	// On distingue  les cas r<>rcut=2.5*sigma
	if (rij <= rm) {
		double f = 48.0 * (eps / rij) * (pow(sig / rij, 12) - ( pow(sig / rij, 6)/2));
		return f;
	}
	else if (rm < rij && rij <= rcut) {
		double x = (rij - rm) / (rcut - rm);
		double f = 48.0 * (eps / rij) * (pow(sig / rij, 12) - (pow(sig / rij, 6) / 2))*(2 * pow(x, 3) - 3 * pow(x, 2) + 1);
		f += -pot(rij)*(6 * pow(x, 2) - 6 * x) / (rcut - rm);
		return f;
	}
	return 0.0;
}

double force2D(double q1[N][2],double F[N][2]) {
	double ep = 0;
	for (int i = 0; i < N; i++) {
		F[i][0] = 0;
		F[i][1] = 0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			double rij[2] = { dist2D(q1[i], q1[j])[0], dist2D(q1[i], q1[j])[1]};

			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2));
			double f = f1D(dist);//Norme de la force

			//Direction ?

			F[i] += -f;
			F[j] += +f;
			ep += pot(dist);
			//cout << ep << endl;
		}
	}
	return ep;
}


double Verlet(double* x0, double* x1, double* v, double* F,double ep){
	double sumv = 0;
	double sumvquad = 0;
	double x2[N];
	double etot = 0;

	for (int i = 0; i < N; i++) {
		x2[i] = 2 * x1[i] - x0[i] + pow(deltaT, 2)*F[i]; //On met � jour les positions
		v[i]= (x2[i] - x0[i]) / (2 * deltaT);//On met � jour les vitesses
		//cout << v[i] << endl;
		sumv += v[i];
		sumvquad += pow(v[i], 2);
		x0[i] = x1[i];
		x1[i] = x2[i];
	}
	/*
	for (int i = 0; i < N; i++) {
		v[i] += -sumv/N; //On fixe la moyenne � 0
	}
	*/
	etot = ep + 0.5*m*sumvquad;
	return etot;
}



int main(){
	InitRandom();

	double q0[N][2];
	double q1[N][2];
	double v[N][2];
	double F[N][2];
	init2D(q0,q1, v);
	force2D(q1, F);
	
	for (int i = 0; i < N; i++) {
		cout << i << ". Position: " << q1[i] << ", vitesse: " << v[i] << ", force: " << F[i] << endl;
	}
	
	//D�marrer l'algorithme
	double t = 0;
	int tour = 0;
	

	ofstream positions;
	positions.open("Positions.txt");
	positions << "Positions des particules \n";
	
	vector< tuple <double,double,double>> Energies;

	//On calcul l'erreur sur l'�nergie totale
	double etot_init = Verlet(x0, x1, v, F, force1D(x1, F));
	double etot_err = 0;

	while (t < tmax) {
		double ep = force1D(x1, F);
		double etot=Verlet(x0, x1, v, F, ep);
		
		etot_err += abs(etot - etot_init);

		if (tour%int(tau/deltaT)==0){
			//On �chantillone tous les tau = 10*deltaT
			// On ins�re les �nergies dans un vector: elles sont dans l'ordre inverse
			Energies.push_back(make_tuple(etot, ep, etot - ep));
			for (int i=0; i < N; i ++) {
				positions << x1[i] <<" ";
			}
			positions << "\n";
		}
		t += deltaT;
		tour++;
	}

	etot_err /= (tmax / deltaT); //On divise par le nombre de valeurs
	cout << "Erreur sur l'energie totale en pourcentage par rapport a l'energie initiale: " << etot_err/etot_init << endl;

	
	
	//Test d�riv�e du potentiel -> force
	//On enregistre les valeurs dans un fichier .txt
	ofstream myfile;
	myfile.open("Derivee.txt");
	double start = 1;
	double step = (rcut - start) / 100;
	myfile << "Distance Pot-prime Force \n";
	for (double i = start; i < rcut; i += step) {
		myfile << to_string(i)+" "+ to_string(-pot_prime(i)) + " " + to_string(f1D(i))+" \n";
	}


	//Exportation de l'�nergie
	ofstream myfile2;
	myfile2.open("Energie.txt");
	myfile2 << "En_tot Ep Ec \n";
	for (vector<tuple<double,double,double>>::iterator it = Energies.begin(); it != Energies.end(); ++it) {
		myfile2 << to_string(get<0>(*it)) + " " + to_string(get<1>(*it)) + " " + to_string(get<2>(*it)) + " \n";
	}
	
	


	return 0;
	}
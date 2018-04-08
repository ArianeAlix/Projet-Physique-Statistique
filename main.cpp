#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
using namespace std;


//Paramètres de la simulation
const int N = 10;//nombre de particules
const double T = 120;//température du système
const double k = 1.3806485279*pow(10, -23);//constante de boltzmann

//Paramètres des particules
const double sigR = 3.405*pow(10, -10);//sigma réelle
const double sig = 1;//sigma réduit
const double a = sig;//taille d'une cellules en sigma
const double epsR = k * T; //vrai epsilon
const double eps = 1;//eps réduit
const double mR = 6.63*pow(10, -26);//masse réelle d'un atome d'argon
const double m = 1;//masse réduite
const double L = N * sig* a;

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.5*sig;

// Constantes liées au temps
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


double dist1D(double xi, double xj) {
	//Doit prendre en compte la périodicité
	//cout << xj - xi << endl;
	double rij = xj - xi - L * round(((xj - xi) / L));
	return rij;
}


void init1D(double* x0,double* x1, double * v,double deltaT) {
	//On choisit un segment de longueur N*a*sigma
	//On place les particules de manière régulière sur le segment
	//Et de façon à ce que la périodicité implique une régularité 
	//des distances entre toutes les particules
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
	double en = 0.0;//energy of the whole system
	double sumv = 0.0;
	double sumvquad = 0.0;
	for (int i = 0; i < N; i++) {
		x1[i] = i * L/N + L/(2.0*N);
		v[i] = distribution(generator);
		sumv += v[i];
		sumvquad += pow(v[i], 2);
	}
	double moyv = sumv / N;
	double moyvquad = sumvquad / N;
	for (int i = 0; i < N; i++) {
		v[i] += -moyv; //On fixe la moyenne à 0
		x0[i] = x1[i] - v[i] * deltaT;
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
	double h = pow(10, -8);
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


double force1D(double* x1,double* F) {
	double ep = 0;
	for (int i = 0; i < N; i++) {
		F[i] = 0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			double rij = dist1D(x1[i], x1[j]);
			double f = f1D(abs(rij));
			F[i] += -f * (rij/abs(rij));
			F[j] += +f * (rij / abs(rij));
			ep += pot(rij);
			//cout << ep << endl;
		}
	}
	return ep;
}


double Verlet(double* x0, double* x1, double* v, double* F,double ep,double deltaT){
	double sumv = 0;
	double sumvquad = 0;
	double x2[N];
	double etot = 0;

	for (int i = 0; i < N; i++) {
		x2[i] = 2 * x1[i] - x0[i] + pow(deltaT, 2)*F[i]; //On met à jour les positions
		v[i]= (x2[i] - x0[i]) / (2 * deltaT);//On met à jour les vitesses
		sumv += v[i];
		sumvquad += pow(v[i], 2);
		x0[i] = x1[i];
		x1[i] = x2[i];
	}
	etot = ep + 0.5*m*sumvquad;
	return etot;
}


//Ensemble de la simulation appelant toutes les fonctions nécessaires
double simulation(double deltaT) { //Renvoie l'erreur
	//Initialisation
	double x0[N];
	double x1[N];
	double v[N];
	double F[N];
	init1D(x0, x1, v, deltaT);
	force1D(x1, F);

	for (int i = 0; i < N; i++) {
		cout << i << ". Position: " << x1[i] << ", vitesse: " << v[i] << ", force: " << F[i] << endl;
	}

	//Pour exportation des positions
	ofstream positions;
	positions.open("Positions.txt");
	positions << "Positions des particules \n";

	//Exportation de l'énergie
	ofstream energies;
	energies.open("Energie.txt");
	energies << "En_tot Ep Ec \n";

	double etot_init = Verlet(x0, x1, v, F, force1D(x1, F), deltaT);
	double etot_err = 0;

	//Démarrer l'algorithme
	double t = 0;
	int tour = 0;

	while (t < tmax) {
		double ep = force1D(x1, F);
		double etot = Verlet(x0, x1, v, F, ep, deltaT);

		etot_err += abs(etot - etot_init);

		if (tour%int(10) == 0) {
			//On échantillonne tous les tau = 10*deltaT
			energies << etot << " " << ep << " "<< etot - ep << "\n";
			for (int i = 0; i < N; i++) {
				positions << x1[i] << " ";
			}
			positions << "\n";
		}
		t += deltaT;
		tour++;
	}

	//On calcul l'erreur sur l'énergie totale
	etot_err /= (tmax / deltaT); //On divise par le nombre de valeurs
	etot_err /= etot_init;

	return etot_err;
}



int main(){
	InitRandom();


	//Test dérivée du potentiel -> force
	//On enregistre les valeurs dans un fichier .txt
	ofstream myfile;
	myfile.open("Derivee.txt");
	double start = 1;
	double step = (rcut - start) / 100;
	myfile << "Distance Pot-prime Force \n";
	for (double i = start; i < rcut; i += step) {
		myfile << to_string(i) + " " + to_string(-pot_prime(i)) + " " + to_string(f1D(i)) + " \n";
	}


	//Lancement de la simulation
	double deltaT = pow(10, -2);//pas de temps de la simulation
	double etot_err=simulation(deltaT);
	
	cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;


	return 0;
	}
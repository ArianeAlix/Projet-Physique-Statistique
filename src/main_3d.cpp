#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
#include<math.h>
using namespace std;


//Paramètres de la simulation
const int N = 125;//->nombre de particules = 125 : 5*5*5
const double T = 120;//température du système
const double kB = 1.3806485279*pow(10, -23);//constante de boltzmann
const double k = 1; //réduit

//Paramètres des particules
const double sigR = 3.405*pow(10, -10);//sigma réelle
const double sig = 1;//sigma réduit
const double a = sig;//taille d'une cellules en sigma
const double epsR = kB * T; //vrai epsilon
const double eps = 1;//eps réduit
const double mR = 6.63*pow(10, -26);//masse réelle d'un atome d'argon
const double m = 1;//masse réduite
const double L = pow(N,1/3) * sig* a; //->côté d'un carré

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.8*sig;

// Constantes liées au temps
const double tmax = 100; //temps de la simulation



double *dist3d(double* qi, double* qj){
	double dist[3];
	//Doit prendre un compte la périodicité

	dist[0] = qj[0] - qi[0] - L * round((qj[0] - qi[0]) / L);
	dist[1] = qj[1] - qi[1] - L * round((qj[1] - qi[1]) / L);
	dist[2] = qj[2] - qi[2] - L * round((qj[2] - qi[2]) / L);
	return dist;
}


void init3d(double q0[N][3] , double q1[N][3], double v[N][3], double deltaT) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
	double en = 0.0;//energy of the whole system
	double sumv[3] = { 0, 0, 0 };
	double sumvquad = 0;

	double nb = pow(N, 1 / 3);//Nb de particules sur un coté

	for (int i = 0; i < N; i++) {
		q1[i][0] = (i % int(nb)) * L/ nb + L/(2.0*nb);
		q1[i][1] = i / int(nb)  * L / nb + L / (2.0*nb);
		q1[i][2] = i % int(pow(nb,2))  * L / nb + L / (2.0*nb);

		v[i][0] = distribution(generator);
		v[i][1] = distribution(generator);
		v[i][2] = distribution(generator);

		sumv[0] += v[i][0];
		sumv[1] += v[i][1];
		sumv[2] += v[i][2];

		sumvquad += pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2);
	}
	double moyv[3] = { sumv[0] / N, sumv[1] / N ,sumv[2] / N };
	double moyvquad = sumvquad / N;

	for (int i = 0; i < N; i++) {
		v[i][0] += -moyv[0]; 
		v[i][1] += -moyv[1]; 
		v[i][2] += -moyv[2];//On fixe la moyenne à 0

		q0[i][0] = q1[i][0] - v[i][0] * deltaT;
		q0[i][1] = q1[i][1] - v[i][1] * deltaT;
		q0[i][2] = q1[i][2] - v[i][2] * deltaT;
	}
	/*
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			double dist[2] = { dist3d(q1[i], q1[j])[0], dist3d(q1[i], q1[j])[1] };
			cout << i << "," << j << ": "<< dist[0] << " , "<< dist[1]<<endl;
		}
	}
	*/
	
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


double force3d(double q1[N][3],double F[N][3]) {
	double ep = 0;
	for (int i = 0; i < N; i++) {
		F[i][0] = 0;
		F[i][1] = 0;
		F[i][2] = 0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			double rij[3] = { dist3d(q1[i], q1[j])[0], dist3d(q1[i], q1[j])[1],dist3d(q1[i], q1[j])[2] };
			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2) + pow(rij[1], 3));
			double f = f1D(dist);//Valeur de la force
			//cout << dist << " ; " << f << endl;
			//Direction ?
			
			F[i][0] += +f * rij[0] / dist;
			//cout << rij[0] << " : " <<dist << " : "<< f <<endl;
			F[i][1] += +f * rij[1] / dist;
			F[i][2] += +f * rij[2] / dist;

			F[j][0] += -f * rij[0] / dist;
			F[j][1] += -f * rij[1] / dist;
			F[j][2] += -f * rij[2] / dist;
			
			
			ep += pot(dist);
			//cout << ep << endl;
		}
		//cout << F[i][0] << " , " << F[i][1] << endl;
	}
	return ep;
}



double Verlet(double q0[N][3], double q1[N][3], double v[N][3], double F[N][3], double ep, double deltaT) {
	double sumv[3] = { 0,0,0 };
	double sumvquad = 0;
	double q2[N][3];
	double etot = 0;

	for (int i = 0; i < N; i++) {
		
		//cout << F[i][0] << " , " << F[i][1] << endl;
		q2[i][0] = 2 * q1[i][0] - q0[i][0] + pow(deltaT, 2)*F[i][0]; //On met à jour les positions
		q2[i][1] = 2 * q1[i][1] - q0[i][1] + pow(deltaT, 2)*F[i][1];
		q2[i][2] = 2 * q1[i][2] - q0[i][2] + pow(deltaT, 2)*F[i][2];

		v[i][0] = (q2[i][0] - q0[i][0]) / (2 * deltaT);//On met à jour les vitesses
		v[i][1] = (q2[i][1] - q0[i][1]) / (2 * deltaT);
		v[i][2] = (q2[i][2] - q0[i][2]) / (2 * deltaT);

		sumv[0] += v[i][0];
		sumv[1] += v[i][1];
		sumv[2] += v[i][2];
		//cout << sumv[0] << " , " << sumv[1] << endl;

		sumvquad += pow(v[i][0], 2) + pow(v[i][1], 2) + pow(v[i][2], 2);

		q0[i][0] = q1[i][0];
		q0[i][1] = q1[i][1];
		q0[i][2] = q1[i][2];

		q1[i][0] = q2[i][0];
		q1[i][1] = q2[i][1];
		q1[i][2] = q2[i][2];
	}
	//cout << q0[0][0] << " , " << q0[0][1] << " : " << F[0][0] << " , " << F[0][1] <<endl;
	//cout << q0[1][0] << " , " << q0[1][1] << " : " << F[1][0] << " , " << F[1][1] << endl;
	etot = ep + 0.5*m*sumvquad;
	return etot;
}


//Ensemble de la simulation appelant toutes les fonctions nécessaires
double simulation(double deltaT) { //Renvoie l'erreur
								   //Initialisation
	double q0[N][3];
	double q1[N][3];
	double v[N][3];
	double F[N][3];
	init3d(q0, q1, v, deltaT);
	force3d(q1, F);

	//Pour afficher les valeurs initiales
	
	for (int i = 0; i < N; i++) {
	cout << i << ". Position: " << q1[i][0] <<"," << q1[i][1] << "," << q1[i][2] << ", vitesse: " << v[i][0] << "," << v[i][1] << "," << v[i][2] << ", force: " << F[i][0] << "," << F[i][1] << "," << F[i][2] << endl;
	}
	

	//Pour exportation des positions
	ofstream positions;
	positions.open("Positions_3d.txt");
	positions << "Positions des particules \n";

	//Exportation de l'énergie
	ofstream energies;
	energies.open("Energie_3d.txt");
	energies << "En_tot Ep Ec \n";

	//Exportation de la température
	ofstream temperatures;
	temperatures.open("Temp_3d.txt");
	temperatures << "Tc Tp \n";


	double etot_init = Verlet(q0, q1, v, F, force3d(q1, F), deltaT);
	double etot_err = 0;

	//Démarrer l'algorithme
	double t = 0;
	int tour = 0;

	while (t < 0.04) {
		double ep = force3d(q1, F);
		double etot = Verlet(q0, q1, v, F, ep, deltaT);

		etot_err += abs(etot - etot_init);

		if (tour%int(10) == 0) {
			//On échantillonne tous les tau = 10*deltaT
			//cout << ep << endl;
			energies << etot << " " << ep << " " << etot - ep << "\n";
			temperatures << 2 * (etot - ep) / (N*k) << " " << 2 * ep / (N*k) << "\n";
			for (int i = 0; i < N; i++) {
				positions << q1[i][0] << " " << q1[i][1] << " ";
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

/*

void erreur_energie() {
	//Lancer plusieurs simulations avec des deltaT différents
	//pour tracer l'erreur en fonction de deltaT^2
	ofstream erreur_en;
	erreur_en.open("Erreur_energie_3d.txt");
	erreur_en << "deltaT etot_err \n";

	double deltaT = pow(10, -5);
	int compteur = 0;
	for (int i = 0; i < 30; i++) {
		cout << "Simulation numero " << ++compteur << "..." << endl;
		deltaT += 5 * pow(10, -4);//pas de temps de la simulation
		erreur_en << pow(deltaT, 2) << " " << simulation(deltaT) << "\n";
	}
	cout << "Fin" << endl;
}

*/

int main() {

	double q0[N][3];
	double q1[N][3];
	double v[N][3];
	double F[N][3];

	
	//Lancement d'une simulation
	double deltaT = pow(10, -3);//pas de temps de la simulation
	double etot_err = simulation(deltaT);

	cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;

	/*
	//Pour pouvoir tracer l'erreur sur l'énergie en fonction du pas de temps
	erreur_energie();

	*/
	return 0;
}
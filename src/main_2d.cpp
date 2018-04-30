#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
#include<math.h>
using namespace std;


//Param�tres de la simulation
const int N = 36;//->nombre de particules = 100 : 10*10
const double T = 120;//temp�rature du syst�me
const double kB = 1.3806485279*pow(10, -23);//constante de boltzmann
const double k = 1; //r�duit

//Param�tres des particules
const double sigR = 3.405*pow(10, -10);//sigma r�elle
const double sig = 1;//sigma r�duit
const double a = sig;//taille d'une cellules en sigma
const double epsR = kB * T; //vrai epsilon
const double eps = 1;//eps r�duit
const double mR = 6.63*pow(10, -26);//masse r�elle d'un atome d'argon
const double m = 1;//masse r�duite
const double L = sqrt(N) * sig* a; //->c�t� d'un carr�

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.8*sig;

// Constantes li�es au temps
const double tmax = 0.04; //temps de la simulation


double* dist2D(double* qi, double* qj){
	double dist[2];
	//Doit prendre un compte la p�riodicit�

	dist[0] = qj[0] - qi[0] - L * round((qj[0] - qi[0]) / L);
	dist[1] = qj[1] - qi[1] - L * round((qj[1] - qi[1]) / L);
	return dist;
}


void init2D(double (&q1)[N][2], double (&p)[N][2], double deltaT) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
	double en = 0.0;//energy of the whole system
	double sump[2] = { 0,0 };
	double sumpquad = 0;
	for (int i = 0; i < N; i++) {
		q1[i][0] = (i % int(sqrt(N))) * L/sqrt(N) + L/(2.0*sqrt(N));
		q1[i][1] = i / int(sqrt(N))  * L/sqrt(N) + L/(2.0*sqrt(N));

		p[i][0] = distribution(generator)*m;
		p[i][1] = distribution(generator)*m;

		sump[0] += p[i][0];
		sump[1] += p[i][1];
		sumpquad += pow(p[i][0], 2) + pow(p[i][1], 2);
	}
	double moyp[2] = { sump[0] / N, sump[1] / N };
	double moypquad = sumpquad / N;

	for (int i = 0; i < N; i++) {
		p[i][0] += -moyp[0]; 
		p[i][1] += -moyp[1]; //On fixe la moyenne � 0

	}
	/*
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			double dist[2] = { dist2D(q1[i], q1[j])[0], dist2D(q1[i], q1[j])[1] };
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
		f += -4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6))*(6 * pow(x, 2) - 6 * x) / (rcut - rm);
		return f;
	}
	return 0.0;
}


double force2D(double(&q1)[N][2],double(&F)[N][2]) {
	double ep = 0;
	for (int i = 0; i < N; i++) {
		F[i][0] = 0;
		F[i][1] = 0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			double rij[2] = { dist2D(q1[i], q1[j])[0], dist2D(q1[i], q1[j])[1]};
			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2));
			double f = f1D(dist);//Valeur de la force
			//cout << dist << " ; " << f << endl;
			//Direction ?
			
			F[i][0] += +f * rij[0] / dist;
			//cout << rij[0] << " : " <<dist << " : "<< f <<endl;
			F[i][1] += +f * rij[1] / dist;
			F[j][0] += -f * rij[0] / dist;
			F[j][1] += -f * rij[1] / dist;
			
			
			ep += pot(dist);
			//cout << ep << endl;
		}
		//cout << F[i][0] << " , " << F[i][1] << endl;
	}
	return ep;
}



double Verlet(double(&q1)[N][2], double(&p)[N][2], double(&F)[N][2], double ep, double deltaT) {
	double sump[2] = { 0,0 };
	double sumpquad = 0;
	double q2[N][2];
	double etot = 0;

	double p_inter[N][2];

	for (int i = 0; i < N; i++) {
		sump[0] += p[i][0];
		sump[1] += p[i][1];
		sumpquad += pow(p[i][0], 2) + pow(p[i][1], 2);

		p_inter[i][0] = p[i][0] + (deltaT / 2)*F[i][0];//quantit� de mvt interm�diaire 
		p_inter[i][1] = p[i][1] + (deltaT / 2)*F[i][1];

		q2[i][0] = q1[i][0] + deltaT * (p_inter[i][0] / m);
		q2[i][1] = q1[i][1] + deltaT * (p_inter[i][1] / m);

	}

	ep = force2D(q2, F); //On recalcule les forces pour les nouvelles positions

	for (int i = 0; i < N; i++) {
		p[i][0] = p_inter[i][0] + (deltaT / 2)*F[i][0];//On met � jour la quantit� de mvt avec les nouvelles forces
		p[i][1] = p_inter[i][1] + (deltaT / 2)*F[i][1];

		q1[i][0] = q2[i][0];//On enregistre les nouvelles positions
		q1[i][1] = q2[i][1];
	}

	etot = ep + 0.5*sumpquad / m;
	return etot;
}


//Ensemble de la simulation appelant toutes les fonctions n�cessaires
double simulation(double deltaT) { //Renvoie l'erreur
								   //Initialisation
	double q1[N][2];
	double p[N][2];
	double F[N][2];
	init2D(q1, p, deltaT);
	force2D(q1, F);

	//Pour afficher les valeurs initiales
	
	for (int i = 0; i < N; i++) {
	cout << i << ". Position: " << q1[i][0] <<"," << q1[i][1] << ", vitesse: " << p[i][0] << "," << p[i][1] << ", force: " << F[i][0] << "," << F[i][1] << endl;
	}
	

	//Pour exportation des positions
	ofstream positions;
	positions.open("Positions_2d.txt");
	positions << "Positions des particules \n";

	//Exportation de l'�nergie
	ofstream energies;
	energies.open("Energie_2d.txt");
	energies << "En_tot Ep Ec \n";

	//Exportation de la temp�rature
	ofstream temperatures;
	temperatures.open("Temp_2d.txt");
	temperatures << "Tc Tp \n";


	double etot_init = Verlet( q1, p, F, force2D(q1, F), deltaT);
	double etot_err = 0;

	//D�marrer l'algorithme
	double t = 0;
	int tour = 0;

	while (t < tmax) {
		double ep = force2D(q1, F);
		double etot = Verlet( q1, p, F, ep, deltaT);

		etot_err += abs(etot - etot_init);

		if (tour%int(10) == 0) {
			//On �chantillonne tous les tau = 10*deltaT
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

	//On calcul l'erreur sur l'�nergie totale
	etot_err /= (tmax / deltaT); //On divise par le nombre de valeurs
	etot_err /= etot_init;

	return etot_err;
}

/*

void erreur_energie() {
	//Lancer plusieurs simulations avec des deltaT diff�rents
	//pour tracer l'erreur en fonction de deltaT^2
	ofstream erreur_en;
	erreur_en.open("Erreur_energie_2d.txt");
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


void perturbation(double(&q1)[N][2], double delt[N][2]) {
	double q2[N][2];
	double q3[N][2];

	for (int i = 0; i < N; i++) {
		//On copie la configuration initiale en 2 nouvelles config, avec la perturbation en delta
		//Positive pour l'une, n�gative pour l'autre
		q2[i][0] = q1[i][0] - delt[i][0];
		q2[i][1] = q1[i][1] - delt[i][1];

		q3[i][0] = q1[i][0] + delt[i][0];
		q3[i][1] = q1[i][1] + delt[i][1];
	}

	// On calcule la diff�rence de potentiel pour chaque particule
	double diff_pot=0;

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			double rij2[2] = { dist2D(q2[i], q2[j])[0], dist2D(q2[i], q2[j])[1] };
			double rij3[2] = { dist2D(q3[i], q3[j])[0], dist2D(q3[i], q3[j])[1] };

			double dist2 = sqrt(pow(rij2[0], 2) + pow(rij2[1], 2));
			double dist3 = sqrt(pow(rij3[0], 2) + pow(rij3[1], 2));

			diff_pot += pot(dist3) - pot(dist2);
		}
	}

	double F1[N][2];
	force2D(q1, F1);
	//On calcule la force qui s'applique pour 
	//en d�duire le gradient * delta sur chaque particule

	double grad_pot_delt=0;

	for (int i = 0; i < N; i++) {
		grad_pot_delt += -(F1[i][0] * delt[i][0] + F1[i][1] * delt[i][1]);
	}

	cout << diff_pot -2* grad_pot_delt << endl;

}

int main() {
	
	//Lancement d'une simulation
	double deltaT = pow(10, -3);//pas de temps de la simulation
	double etot_err = simulation(deltaT);

	cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;

	/*
	//Pour pouvoir tracer l'erreur sur l'�nergie en fonction du pas de temps
	erreur_energie();

	*/

	double q1[N][2];
	double v[N][2];
	double F[N][2];
	init2D(q1, v, deltaT);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, 1);
	double delt[N][2];
	double delta[N][2];
	for (int i = 0; i < N; i++) {
		delt[i][0] = distribution(generator)*pow(10, -6);
		delt[i][1] = distribution(generator)*pow(10, -6);
		delta[i][0] = distribution(generator)*pow(10, -1);
		delta[i][1] = distribution(generator)*pow(10, -1);
	}

	for (int i = 0; i < N; i++) {
		//On copie la configuration initiale en 2 nouvelles config, avec la perturbation en delta
		//Positive pour l'une, n�gative pour l'autre
		q1[i][0] = q1[i][0] - delta[i][0];
		q1[i][1] = q1[i][1] - delta[i][1];
	}

	perturbation(q1, delt);



	return 0;
}
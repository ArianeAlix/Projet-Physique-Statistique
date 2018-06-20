#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
#include<math.h>
using namespace std;


//Paramètres de la simulation
const int N = 36;//->nombre de particules = 100 : 10*10
const double T = 120;//température du système
const double kB = 1.3806485279*pow(10, -23);//constante de boltzmann
const double k = 1; //réduit


//Paramètres des particules
const double sigR = 3.405*pow(10, -10);//sigma réelle
const double sig = 1.0;//sigma réduit
const double a = 1*sig;//taille d'une cellules en sigma
const double epsR = kB * T; //vrai epsilon
const double eps = 1;//eps réduit
const double mR = 6.63*pow(10, -26);//masse réelle d'un atome d'argon
const double m = 1;//masse réduite
const double L = sqrt(N) * a; //->côté d'un carré

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.8*sig;

// Constantes liées au temps
const double tmax = 20; //temps de la simulation

//Pour l'isocinétique
const double T_cible = 1.2;
const double Ec_cible = 2 / 2 * N*eps*T_cible;


void dist2D(double* qi, double* qj, double(&dist)[2]){
	//Doit prendre un compte la périodicité
	dist[0] = qj[0] - qi[0] - L * round((qj[0] - qi[0]) / L);
	dist[1] = qj[1] - qi[1] - L * round((qj[1] - qi[1]) / L);
}


void init2D(double (&q1)[N][2], double (&p)[N][2]) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator(0); //Si l'on veut tracer des paramètres/valeur pour une même simulation
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
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
		p[i][1] += -moyp[1]; //On fixe la moyenne à 0

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

double derivee2(double rij) {// en 1D correspond à la dérivée seconde
							 // On distingue  les cas r<>rcut=2.5*sigma
	if (rij <= rm) {
		double f = 48.0 * (eps / pow(rij, 2)) * (13 * pow(sig / rij, 12) - (7 / 2)*(pow(sig / rij, 6)));
		return f;
	}
	else if (rm < rij && rij <= rcut) {
		double x = (rij - rm) / (rcut - rm);
		double f = 48.0 * (eps / pow(rij, 2)) * (13 * pow(sig / rij, 12) - (7 / 2)*(pow(sig / rij, 6)))*(2 * pow(x, 3) - 3 * pow(x, 2) + 1);
		f += 2 * (-48.0 * (eps / rij) * (pow(sig / rij, 12) - (pow(sig / rij, 6) / 2))*(6 * pow(x, 2) - 6 * x) / (rcut - rm));
		f += 4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6))*(12 * x - 6) / pow((rcut - rm), 2);
		return f;
	}
	return 0.0;
}


//return the total potential of the configuration q1
double potential_tot(double(&q1)[N][2]) {
	double potential = 0;
	double rij[2];
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			dist2D(q1[i], q1[j],rij);
			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2));
			potential += pot(dist);
		}
	}
	return potential;
}


//return the laplacian
double laplacian(double(&q1)[N][2]) {
	double laplace = 0;
	double dV[N];
	double rij[2];
	for (int i = 0; i < N; i++) {
		dV[i] = 0;
	}

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			dist2D(q1[i], q1[j], rij);
			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2));
            dV[i] += derivee2(dist) + f1D(dist)*1/dist;
            dV[j] += derivee2(dist) + f1D(dist)*1 / dist;
		}
	}
	for (int i = 0; i < N; i++) {
		laplace += dV[i];
	}

	return laplace;
}


//F equals to the vector -grad
double force2D(double(&q1)[N][2],double(&F)[N][2], double& moygradquad, double& moylaplace, double & press) {
	double ep = 0;
	press = 0;
	moygradquad = 0;
	moylaplace = 0;

	for (int i = 0; i < N; i++) {
		F[i][0] = 0;
		F[i][1] = 0;
	}
	
	double rij[2];

	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			dist2D(q1[i], q1[j], rij);
			double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2));
			double f = f1D(dist);//Valeur de la force
			//cout << dist << " ; " << f << endl;
	
			F[i][0] += -f * rij[0] / dist;
			F[i][1] += -f * rij[1] / dist;
			F[j][0] += +f * rij[0] / dist;
			F[j][1] += +f * rij[1] / dist;
			
			ep += pot(dist);
			//cout << ep << endl;

			//Calcul pour la pression
			press += dist*f;
			
		}
		//cout << F[i][0] << " , " << F[i][1] << endl;
	}
	for (int i = 0; i < N; i++) {
		moygradquad += pow(F[i][0], 2)+ pow(F[i][1], 2);
	}


	press *= 0.5;
	moygradquad /= N;
	moylaplace =laplacian(q1)/ N;
	return ep;
}



double Verlet(double(&q1)[N][2], double(&p)[N][2], double& ep, double& ec, double& press, double& moygradquad, double& moylaplace, double(&F)[N][2], double deltaT) {
	double sump[2] = { 0,0 };
	double sumpquad = 0;
	double q2[N][2];

	double p_inter[N][2];


	//Facteur à utiliser pour l'isocinétique
	//double lambda = potential_tot(q1) / (2 * Ec_cible);
	

	for (int i = 0; i < N; i++) {
		p_inter[i][0] = p[i][0]  + (deltaT / 2)*F[i][0];//quantité de mvt intermédiaire 
		p_inter[i][1] = p[i][1] + (deltaT / 2)*F[i][1];

		q2[i][0] = q1[i][0] + deltaT * (p_inter[i][0] / m);
		q2[i][1] = q1[i][1] + deltaT * (p_inter[i][1] / m);

	}

	ep = force2D(q2, F, moygradquad, moylaplace, press); //On recalcule les forces pour les nouvelles positions

	for (int i = 0; i < N; i++) {
		p[i][0] = p_inter[i][0] + (deltaT / 2)*F[i][0];//On met à jour la quantité de mvt avec les nouvelles forces
		p[i][1] = p_inter[i][1] + (deltaT / 2)*F[i][1];

		sump[0] += p[i][0];
		sump[1] += p[i][1];
		sumpquad += pow(p[i][0], 2) + pow(p[i][1], 2);

		q1[i][0] = q2[i][0];//On enregistre les nouvelles positions
		q1[i][1] = q2[i][1];
	}

	ec = 0.5*sumpquad / m;

	
	double lambda = sqrt(Ec_cible / ec);
	for (int i = 0; i < N; i++) {
		p[i][0] = p[i][0] * lambda;
		p[i][1] = p[i][1] * lambda;
	}
	
	ec = 0;
	for (int i = 0; i < N; i++) {
		ec += 0.5*(pow(p[i][0], 2) + pow(p[i][1], 2))/m;
	}


	return ep+ec;
}



//Ensemble de la simulation appelant toutes les fonctions nécessaires
double simulation(double deltaT, double& emin, double& emax) { //Renvoie l'erreur
															   //Initialisation
	double q1[N][2];
	double p[N][2];
	double F[N][2];
    init2D(q1, p);
	double moygradquad;
	double moylaplace;
	double ec;
	double press;
	double pression_sum = 0;

	//Pour afficher les valeurs initiales
	/*
	for (int i = 0; i < N; i++) {
		cout << i << ". Position: " << q1[i][0] << "," << q1[i][1] << ", vitesse: " << p[i][0] << "," << p[i][1] << ", force: " << F[i][0] << "," << F[i][1] << endl;
	}
	*/
	double ep = force2D(q1, F, moygradquad, moylaplace, press);


	//Pour exportation des positions
	ofstream positions;
	positions.open("Positions_2d.txt");
	positions << "Positions des particules \n";

	//Exportation de l'énergie
	ofstream energies;
	energies.open("Energie_2d.txt");
	energies << "En_tot Ep Ec \n";

	//Exportation de la température
	ofstream temperatures;
	temperatures.open("Temp_2d.txt");
	temperatures << "Tc Tp \n";

	//Exportation de la pression
	ofstream pression;
	pression.open("Press_2d.txt", std::ios_base::app | std::ios_base::out);
	//pression << "Densite P \n";

	//Exportation de la pression
	ofstream Epmoy;
	Epmoy.open("Ep_2d.txt", std::ios_base::app | std::ios_base::out);
	//pression << "Densite P \n";
	double epmoy = 0;

	//Histogramme des vitesses pour une particule
	ofstream Histogramme;
	Histogramme.open("Histo_2d.txt", std::ios_base::app | std::ios_base::out);



	double etot_init = Verlet(q1, p, ep, ec ,press,moygradquad, moylaplace, F, deltaT);
	emin = etot_init;
	emax = etot_init;
	double etot_err = 0;

	//Démarrer l'algorithme
	double t = 0;
	int tour = 0;

	

	while (t < tmax) {
		//double ep = force1D(x1, F);
		double etot = Verlet(q1, p, ep, ec,press,  moygradquad, moylaplace, F, deltaT);

		if (etot > emax) {
			emax = etot;
		}
		if (etot < emin) {
			emin = etot;
		}

		etot_err += abs(etot - etot_init);
		pression_sum += press; //On ajoute la pression courante 
		epmoy += ep;
		Histogramme << p[0][0] / m << " " << p[0][1] / m << "\n";

		if (tour%int(10) == 0) {
			//On échantillonne tous les tau = 10*deltaT
			//cout << ep << endl;
			energies << etot << " " << ep << " " << ec << "\n";
			temperatures << 2 * ec / (N*k*2) << " " << moygradquad / moylaplace << "\n";
			

			for (int i = 0; i < N; i++) {
				positions << q1[i][0] << " " << q1[i][1] << " ";
			}
			positions << "\n";
		}
		t += deltaT;
		tour++;

	}

	//On calcule l'erreur sur l'énergie totale
	etot_err /= (tmax / deltaT); //On divise par le nombre de valeurs
	etot_err /= etot_init;

	pression << N / pow(L, 2) << " " << (1 / (2 * pow(L, 2)))*(2*N*k*T_cible +     pression_sum / (tmax/deltaT) )<< "\n";
	Epmoy << N / pow(L, 2) << " " << epmoy / (tmax / deltaT) << "\n";

	return etot_err;
}

void erreur_energie() {
	//Lancer plusieurs simulations avec des deltaT différents
	//pour tracer l'erreur en fonction de deltaT^2
	ofstream erreur_en;
	erreur_en.open("Erreur_energie_2d.txt");
	erreur_en << "deltaT etot_err \n";


	double deltaT = pow(10, -4);
	int compteur = 0;


	for (int i = 0; i < 10; i++) {
		double emin, emax;
		cout << "Simulation numero " << ++compteur << "..." << endl;

		simulation(deltaT, emin, emax);

		deltaT += pow(10, -3);//pas de temps de la simulation
		erreur_en << pow(deltaT, 2) << " " << (emax - emin)/N << "\n";
	}
	cout << "Fin" << endl;
}



void perturbation(double(&q1)[N][2], double delt[N][2]) {
	double q2[N][2];
	double q3[N][2];
	double moygradquad;
	double moylaplace;
	double press;

	for (int i = 0; i < N; i++) {
		//On copie la configuration initiale en 2 nouvelles config, avec la perturbation en delta
		//Positive pour l'une, négative pour l'autre
		q2[i][0] = q1[i][0] - delt[i][0];
		q2[i][1] = q1[i][1] - delt[i][1];

		q3[i][0] = q1[i][0] + delt[i][0];
		q3[i][1] = q1[i][1] + delt[i][1];
	}

	// On calcule la différence de potentiel pour chaque particule
	double diff_pot = potential_tot(q3) - potential_tot(q2);

	double F1[N][2];
	force2D(q1, F1,moygradquad,moylaplace, press);
	//On calcule la force qui s'applique pour 
	//en déduire le gradient * delta sur chaque particule

	double grad_pot_delt=0;

	for (int i = 0; i < N; i++) {
		grad_pot_delt += (-F1[i][0] * delt[i][0] - F1[i][1] * delt[i][1]);
	}

	cout << diff_pot << " , " << 2* grad_pot_delt << endl;

}

int main() {
	
	double emin, emax;
	//Lancement d'une simulation
	double deltaT = pow(10, -4);//pas de temps de la simulation
	double etot_err = simulation(deltaT,emin,emax);

	//cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;
	
	//Pour pouvoir tracer l'erreur sur l'énergie en fonction du pas de temps
	//erreur_energie(); //Déjà fait une fois

	
	//// Validation du gradient en étudiant la différence de potentiel entre deux systèmes perturbés
	/*
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
		//Positive pour l'une, négative pour l'autre
		q1[i][0] = q1[i][0] - delta[i][0];
		q1[i][1] = q1[i][1] - delta[i][1];
	}

	perturbation(q1, delt);
	*/

	return 0;
}

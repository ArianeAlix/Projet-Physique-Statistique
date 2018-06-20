#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
using namespace std;


//Paramètres de la simulation
const int N = 100;//nombre de particules
const double T = 1;
const double TR = 120;//température réelle du système
const double kB = 1.3806485279*pow(10, -23);//constante de boltzmann
const double k = 1; //réduit

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
const double rm = 1.8*sig;

// Constantes liées au temps
const double tmax = 10; //temps de la simulation


double dist1D(double xi, double xj) {
	//Doit prendre en compte la périodicité
	//cout << xj - xi << endl;
	double rij = xj - xi - L * round(((xj - xi) / L));
	return rij;
}


void init1D(double* x1, double * p) {
	//On choisit un segment de longueur N*a*sigma
	//On place les particules de manière régulière sur le segment
	//Et de façon à ce que la périodicité implique une régularité 
	//des distances entre toutes les particules
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//default_random_engine generator(0); //Si l'on veut tracer des paramètres/valeur pour une même simulation
	default_random_engine generator(seed); 
	normal_distribution<double> distribution(0, eps/m);
	double sump = 0.0;
	double sumpquad = 0.0;

	for (int i = 0; i < N; i++) {
		x1[i] = i * L/N + L/(2.0*N);
		p[i] = distribution(generator)*m;
		sump += p[i];
		sumpquad += pow(p[i], 2);
	}
	double moyp = sump / N;
	double moypquad = sumpquad / N;
	for (int i = 0; i < N; i++) {
		p[i] += -moyp; //On fixe la moyenne à 0
	}

	/*
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			double dist = dist1D(x1[i], x1[j]);
			cout << i << "," << j << ": "<< dist <<endl;
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

double pot_real(double rij) {
	double p = 4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6));
	return p;
}


double pot_prime(double rij) {
	double h = pow(10, -8);
	double pp = (pot(rij + h) - pot(rij-h)) / (2*h);
	return pp;
}

double pot_prime_prime(double rij) {
	double h = pow(10, -5);
	double pp = (pot_prime(rij + h) - pot_prime(rij - h)) / (2 * h);
	return pp;
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
		double f = 48.0 * (eps / pow(rij,2)) * (13*pow(sig / rij, 12) - (7/2)*(pow(sig / rij, 6)));
		return f;
	}
	else if (rm < rij && rij <= rcut) {
		double x = (rij - rm) / (rcut - rm);
		double f = 48.0 * (eps / pow(rij, 2)) * (13 * pow(sig / rij, 12) - (7 / 2)*(pow(sig / rij, 6)))*(2 * pow(x, 3) - 3 * pow(x, 2) + 1);
		f += 2*(-48.0 * (eps / rij) * (pow(sig / rij, 12) - (pow(sig / rij, 6) / 2))*(6 * pow(x, 2) - 6 * x) / (rcut - rm));
		f += 4 * eps * (pow(sig / rij, 12) - pow(sig / rij, 6))*(12 * x-6) / pow((rcut - rm), 2);
		return f;
	}
	return 0.0;
}



double force1D(double* x1,double* F,double& moygradquad, double& moylaplace) {
	double ep = 0;
	moygradquad = 0;
	moylaplace = 0;
	double laplace[N];
	for (int i = 0; i < N; i++) {
		F[i] = 0;
		laplace[i] = 0;
	}

	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			double rij = dist1D(x1[i], x1[j]);
			double f = f1D(abs(rij));
			//cout << abs(rij) << " , " <<f <<" , " << derivee2(abs(rij)) << endl;
			laplace[i] += derivee2(abs(rij));
			laplace[j] += derivee2(abs(rij));
			//cout << f << endl;
			F[i] += -f * (rij/abs(rij));
			F[j] += +f * (rij / abs(rij));
			ep += pot(abs(rij));
			//cout << ep << endl;
		}
	}
	for (int i = 0; i < N; i++) {
		moygradquad += pow(F[i], 2);
		//cout << moylaplace << endl;
		moylaplace += laplace[i];
	}

	moygradquad /=N ;
	moylaplace /= N ;
	return ep;
}


double Verlet(double* x1, double* p, double& ep, double& moygradquad, double& moylaplace, double* F, double deltaT){
	double sump = 0;
	double sumpquad = 0;
	double x2[N];
	double etot = 0;
	double p_inter[N];

	for (int i = 0; i < N; i++) {
		sump += p[i];
		sumpquad += pow(p[i], 2);
		
		p_inter[i] = p[i] + (deltaT / 2)*F[i];//quantité de mvt intermédiaire 

		x2[i] = x1[i] + deltaT*(p_inter[i]/m);

	}

	ep = force1D(x2, F, moygradquad, moylaplace); //On recalcule les forces pour les nouvelles positions

	for (int i = 0; i < N; i++) {
		p[i] = p_inter[i]+(deltaT / 2)*F[i];//On met à jour la quantité de mvt avec les nouvelles forces

		x1[i] = x2[i];//On enregistre les nouvelles positions
	}

	etot = ep + 0.5*sumpquad/m;
	return etot;
}


//Ensemble de la simulation appelant toutes les fonctions nécessaires
double simulation(double deltaT,double& emin, double& emax) { //Renvoie l'erreur
	//Initialisation
	double x1[N];
	double p[N];
	double F[N];
    init1D( x1, p);
	double moygradquad;
	double moylaplace;

	double ep= force1D(x1, F,moygradquad,moylaplace);


	//Pour afficher les valeurs initiales
	/*
	for (int i = 0; i < N; i++) {
		cout << i << ". Position: " << x1[i] << ", vitesse: " << v[i] << ", force: " << F[i] << endl;
	}
	*/
	

	//Pour exportation des positions
	ofstream positions;
	positions.open("Positions.txt");
	positions << "Positions des particules \n";

	//Exportation de l'énergie
	ofstream energies;
	energies.open("Energie.txt");
	energies << "En_tot Ep Ec \n";

	//Exportation de la température
	ofstream temperatures;
	temperatures.open("Temp.txt");
	temperatures << "Tc Tp \n";

	double etot_init = Verlet( x1, p, ep, moygradquad, moylaplace, F,deltaT);
	emin = etot_init;
	emax = etot_init;
	double etot_err = 0;
	
	//Démarrer l'algorithme
	double t = 0;
	int tour = 0;

	while (t < tmax) {
		//double ep = force1D(x1, F);
		double etot = Verlet( x1, p, ep, moygradquad, moylaplace, F, deltaT);

		if (etot > emax) {
			emax = etot;
		}
		if (etot < emin) {
			emin = etot;
		}

		etot_err += abs(etot - etot_init);

		if (tour%int(10) == 0) {
			//On échantillonne tous les tau = 10*deltaT
			energies << etot << " " << ep << " "<< etot - ep << "\n";
			temperatures << 2*(etot - ep)/(N*k) << " " << moygradquad/moylaplace << "\n";
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



void erreur_energie() {
	//Lancer plusieurs simulations avec des deltaT différents
	//pour tracer l'erreur en fonction de deltaT^2
	ofstream erreur_en;
	erreur_en.open("Erreur_energie.txt");
	erreur_en << "deltaT etot_err \n";

	double deltaT= pow(10, -4);
	int compteur = 0;
	

	for (int i = 0; i < 10; i ++) {
		double emin, emax;
		cout << "Simulation numero " << ++compteur << "..." << endl;

		simulation(deltaT, emin, emax);

		deltaT += pow(10, -3);//pas de temps de la simulation
		erreur_en << pow(deltaT,2) << " " << (emax-emin) << "\n"; 
	}
	cout << "Fin" << endl;
}



int main(){

	//Test dérivée du potentiel -> force
	//On enregistre les valeurs dans un fichier .txt
	ofstream myfile;
	ofstream myfile2;
	myfile.open("Derivee.txt");
	myfile2.open("Derivee2.txt");
	double start = 1;
	double step = (rcut - start) / 100;
	double err_pot = 0;
	double err_pot2 = 0;
	int iter = 0; //Nb d'itérations pour calculer l'erreur

	myfile << "Distance Pot-prime Force \n";
	myfile2 << "Distance Pot-prime-prime Laplacien \n";

	for (double i = start; i < rcut; i += step) {
		myfile << to_string(i) + " " + to_string(-pot_prime(i)) + " " + to_string(f1D(i)) + " \n";
		err_pot += (abs(f1D(i) + pot_prime(i)))/ pot_prime(i); //erreur relative

		myfile2 << to_string(i) + " " + to_string(pot_prime_prime(i)) + " " + to_string(derivee2(i)) + " \n";
		err_pot2 += (abs(derivee2(i) - pot_prime_prime(i))) / pot_prime_prime(i); //erreur relative
		iter += 1;
	}
	err_pot /= iter;//On divise l'erreur totale par le nombre d'itérations
	err_pot2 /= iter;
	cout << "Erreur relative de la force par rapport a la derivee du potentiel: " << err_pot * 100 << "%" << endl;
	cout << "Erreur relative de la derivee seconde par rapport a la difference finie seconde du potentiel: " << err_pot2 * 100 << "%" << endl;



	//Lancement d'une simulation
	double deltaT = pow(10, -3);//pas de temps de la simulation
	double emin, emax;
	double etot_err=simulation(deltaT,emin,emax);
	
	cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;

	/*
	//Pour pouvoir tracer l'erreur sur l'énergie en fonction du pas de temps
	erreur_energie();
	*/

	return 0;
	}

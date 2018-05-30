#include <iostream>
#include <fstream>
#include<random>
#include<chrono>
#include<string>
#include<tuple>
#include<math.h>
using namespace std;


//Paramètres de la simulation
const int N = 64;//->nombre de particules = 64: 4*4*4
const double T = 120;//température du système
const double kB = 1.3806485279*pow(10, -23);//constante de boltzmann
const double k = 1; //réduit

//Paramètres des particules
const double sigR = 3.405*pow(10, -10);//sigma réelle
const double sig = 1;//sigma réduit
const double a = 1*sig;//taille d'une cellules en sigma
const double epsR = kB * T; //vrai epsilon
const double eps = 1;//eps réduit
const double mR = 6.63*pow(10, -26);//masse réelle d'un atome d'argon
const double m = 1;//masse réduite
const double L = pow(N,1.0/3) * sig* a; //->côté d'un carré

//Dimensionnement du potentiel
const double rcut = 2.5*sig;
const double rm = 1.8*sig;

// Constantes liées au temps
const double tmax = 10; //temps de la simulation



double* dist3D(double* qi, double* qj){
    double* dist;
    dist=(double*)malloc(3);
	//Doit prendre un compte la périodicité

	dist[0] = qj[0] - qi[0] - L * round((qj[0] - qi[0]) / L);
	dist[1] = qj[1] - qi[1] - L * round((qj[1] - qi[1]) / L);
	dist[2] = qj[2] - qi[2] - L * round((qj[2] - qi[2]) / L);
	return dist;
}



void init3D(double (&q1)[N][3], double (&p)[N][3]) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //default_random_engine generator(0); //Si l'on veut tracer des paramètres/valeur pour une même simulation
	default_random_engine generator(seed);
	normal_distribution<double> distribution(0, eps/m);
	double en = 0.0;//energy of the whole system
    double sump[3] = { 0, 0, 0 };
    double sumpquad = 0;

    double nb = pow(N, 1.0 / 3);//Nb de particules sur un coté

	for (int i = 0; i < N; i++) {
        //cout << nb << " " << int(round(nb)) << " " << (i % int(round(nb)))<<endl;
        q1[i][0] = (i % int(round(nb))) *(L/ nb) + L/(2.0*nb);
        q1[i][1] = ((i / int(round(nb))) % int(round(nb))) * (L / nb) + L / (2.0*nb);
        q1[i][2] = (i / int(round(pow(nb,2)))) * (L / nb) + L / (2.0*nb);

        p[i][0] = distribution(generator)*m;
        p[i][1] = distribution(generator)*m;
        p[i][2] = distribution(generator)*m;

        sump[0] += p[i][0];
        sump[1] += p[i][1];
        sump[2] += p[i][2];

        sumpquad += pow(p[i][0], 2) + pow(p[i][1], 2) + pow(p[i][2], 2);
	}
    double moyp[3] = { sump[0] / N, sump[1] / N ,sump[2] / N };
    double moypquad = sumpquad / N;

	for (int i = 0; i < N; i++) {
        p[i][0] += -moyp[0];
        p[i][1] += -moyp[1];
        p[i][2] += -moyp[2];//On fixe la moyenne à 0
	}
    /*
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist[3] = { dist3D(q1[i], q1[j])[0], dist3D(q1[i], q1[j])[1] , dist3D(q1[i], q1[j])[2] };
            cout << i << "," << j << ": "<< dist[0] << " , "<< dist[1]<< " , "<< dist[2]<< endl;
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
double potential_tot(double(&q1)[N][3]) {
    double potential = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double rij[3] = { dist3D(q1[i], q1[j])[0], dist3D(q1[i], q1[j])[1] , dist3D(q1[i], q1[j])[2] };
            double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2)+pow(rij[2], 2));
            potential += pot(dist);
        }
    }
    return potential;
}



//return the laplacian
double laplacian(double(&q1)[N][3]) {
    double laplace = 0;
    double dV[N];

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double rij[3] = { dist3D(q1[i], q1[j])[0], dist3D(q1[i], q1[j])[1] , dist3D(q1[i], q1[j])[2] };
            double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2)+pow(rij[2], 2));
            dV[i] += derivee2(dist) + f1D(dist)* 2/dist;
            dV[j] += derivee2(dist) + f1D(dist)* 2/ dist;
        }
    }
    for (int i = 0; i < N; i++) {
        laplace += dV[i];
    }

    return laplace;
}

//F equals to the vector -grad
double force3D(double(&q1)[N][3],double(&F)[N][3], double& moygradquad, double& moylaplace) {
    double ep = 0;
    moygradquad = 0;
    moylaplace = 0;
    double laplace[N];
    for (int i = 0; i < N; i++) {
        F[i][0] = 0;
        F[i][1] = 0;
        F[i][2] = 0;
        laplace[i] = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            double rij[3] = { dist3D(q1[i], q1[j])[0], dist3D(q1[i], q1[j])[1] , dist3D(q1[i], q1[j])[2] };
            double dist = sqrt(pow(rij[0], 2) + pow(rij[1], 2)+pow(rij[2], 2));
            double f = f1D(dist);//Valeur de la force
            //cout << dist << " ; " << f << endl;
            //Direction ?
            laplace[i] += derivee2(dist);
            laplace[j] += derivee2(dist);

            F[i][0] += -f * rij[0] / dist;
            F[i][1] += -f * rij[1] / dist;
            F[i][2] += -f * rij[2] / dist;

            F[j][0] += +f * rij[0] / dist;
            F[j][1] += +f * rij[1] / dist;
            F[j][2] += +f * rij[2] / dist;

            ep += pot(dist);
            //cout << ep << endl;
        }
        //cout << F[i][0] << " , " << F[i][1] << endl;
    }
    for (int i = 0; i < N; i++) {
        moygradquad += pow(F[i][0], 2)+ pow(F[i][1], 2)+ pow(F[i][2], 2);
        //cout << moylaplace << endl;
        moylaplace += laplace[i];
    }

    moygradquad /= N;
    moylaplace /= N;
    return ep;
}



double Verlet(double(&q1)[N][3], double(&p)[N][3], double& ep, double& moygradquad, double& moylaplace, double(&F)[N][3], double deltaT) {
    double sump[3] = { 0,0,0 };
    double sumpquad = 0;
    double q2[N][3];
    double etot = 0;

    double p_inter[N][3];

    for (int i = 0; i < N; i++) {
        sump[0] += p[i][0];
        sump[1] += p[i][1];
        sump[2] += p[i][2];
        sumpquad += pow(p[i][0], 2) + pow(p[i][1], 2)+ pow(p[i][2], 2);

        p_inter[i][0] = p[i][0] + (deltaT / 2)*F[i][0];//quantité de mvt intermédiaire
        p_inter[i][1] = p[i][1] + (deltaT / 2)*F[i][1];
        p_inter[i][2] = p[i][2] + (deltaT / 2)*F[i][2];

        q2[i][0] = q1[i][0] + deltaT * (p_inter[i][0] / m);
        q2[i][1] = q1[i][1] + deltaT * (p_inter[i][1] / m);
        q2[i][2] = q1[i][2] + deltaT * (p_inter[i][2] / m);

    }

    ep = force3D(q2, F, moygradquad, moylaplace); //On recalcule les forces pour les nouvelles positions

    for (int i = 0; i < N; i++) {
        p[i][0] = p_inter[i][0] + (deltaT / 2)*F[i][0];//On met à jour la quantité de mvt avec les nouvelles forces
        p[i][1] = p_inter[i][1] + (deltaT / 2)*F[i][1];
        p[i][2] = p_inter[i][2] + (deltaT / 2)*F[i][2];

        q1[i][0] = q2[i][0];//On enregistre les nouvelles positions
        q1[i][1] = q2[i][1];
        q1[i][2] = q2[i][2];
    }

    etot = ep + 0.5*sumpquad / m;
    return etot;
}


//Ensemble de la simulation appelant toutes les fonctions nécessaires
double simulation(double deltaT, double& emin, double& emax) { //Renvoie l'erreur
                                                               //Initialisation
    double q1[N][3];
    double p[N][3];
    double F[N][3];
    init3D(q1, p);
    double moygradquad;
    double moylaplace;

    //Pour afficher les valeurs initiales
    /*
    for (int i = 0; i < 10; i++) {
        cout << i << ". Position: " << q1[i][0] << "," << q1[i][1]<< "," << q1[i][2] << ", \n vitesse: " << p[i][0] << "," << p[i][1] << "," << p[i][2] << ", \n force: " << F[i][0] << "," << F[i][1]  << "," << F[i][2]<< endl;
    }
    */

    double ep = force3D(q1, F, moygradquad, moylaplace);


    //Pour exportation des positions
    ofstream positions;
    positions.open("Positions_3d.xyz");

    //Exportation de l'énergie
    ofstream energies;
    energies.open("Energie_3d.txt");
    energies << "En_tot Ep Ec \n";

    //Exportation de la température
    ofstream temperatures;
    temperatures.open("Temp_3d.txt");
    temperatures << "Tc Tp \n";


    double etot_init = Verlet(q1, p, ep, moygradquad, moylaplace, F, deltaT);
    emin = etot_init;
    emax = etot_init;
    double etot_err = 0;

    //Démarrer l'algorithme
    double t = 0;
    int tour = 0;

    while (t < tmax) {
        //double ep = force1D(x1, F);
        double etot = Verlet(q1, p, ep, moygradquad, moylaplace, F, deltaT);

        if (etot > emax) {
            emax = etot;
        }
        if (etot < emin) {
            emin = etot;
        }

        etot_err += abs(etot - etot_init);

        if (tour%int(10) == 0) {
            //On échantillonne tous les tau = 10*deltaT
            //cout << ep << endl;
            energies << etot << " " << ep << " " << etot - ep << "\n";
            temperatures << 2 * (etot - ep) / (N*k*3) << " " << moygradquad / moylaplace << "\n";
            //On doit respecter le format nécessaire pour xmakemol
            positions << to_string(N) << "\n \n";
            for (int i = 0; i < N; i++) {
                positions << "Ar " << q1[i][0] << " " << q1[i][1] << " "  << q1[i][2] << "\n";
            }
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
    double emin, emax;
    //Lancement d'une simulation
    double deltaT = pow(10, -3);//pas de temps de la simulation
    double etot_err = simulation(deltaT,emin,emax);

	cout << "Erreur sur l'energie totale par rapport a l'energie initiale: " << etot_err * 100 << "%" << endl;

	/*
	//Pour pouvoir tracer l'erreur sur l'énergie en fonction du pas de temps
	erreur_energie();

	*/
	return 0;
}

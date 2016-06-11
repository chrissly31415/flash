/*
 * Molecule.cpp
 *
 *  Created on: Jan 9, 2011
 *      Author: loschen
 */

#include <iostream>
#include <math.h>
#include "Molecule.h"
#include "ScfCalc.h"

//Constructor
Molecule::Molecule(int lN, double lchi_s, double lchi_12, double lchi_O,
		double lfi_b) {
	N = lN;
	chi_s = lchi_s;
	chi_12 = lchi_12;
	chi_O = lchi_O;
	fi_b = lfi_b;
	//initialization and declaration of arrays
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = fi_b;
		//fi_z[z] = 0.5;
		fi_zfree[z] = fi_b;
		fi_zads[z] = 0;
		//fi_total[z] = fi_1b + fi_2b;
		//energy parameters
		u_adsorption[z] = 0;
		// its the wall
		//first layer
		alpha[z] = 0.0;
		if (z == 1) {
			u_adsorption[z] = -chi_s;
			//alpha[z] = 0.31;
		}

		if (z == 0) {
			fi_z[z] = 0.0;
			fi_zfree[z] = 0.0;
		}
		alpha_new[z] = 0;
		u_mixing[z] = 0;
		u_total[z] = 0;
		//alpha[z]=-u_mixing[z];
	}

	//initialisation and declaration of "matrices"
	fi = new double[zmax * N];
	fi_free = new double[zmax * N];
	G = new double[zmax * N];
	G_N = new double[zmax * N];
	G_free = new double[zmax * N];
	for (int i = 0; i < zmax * N; i++) {
		//first entry z=0, s=1 and the last N-1 entries are never used
		fi[i] = fi_b / (double) N;
		fi_free[i] = fi_b / (double) N;
		G[i] = 0.0;
		G_N[i] = 0.0;
		G_free[i] = 0.0;
	}

}

// transformation of indices
//   	z1N1 z1N2 z1N3 	     z1N1 z2N1 z3N1 ...
//	 	z2N1 z2N2   ...  ->
//		z3N1 Z3N2   ...
void Molecule::calcG(int k, double temp) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = 1; s <= N; s++) {
			int pos = 0;
			pos = s + z * N - N;
			//initial conditions
			if (s == 1 && k == 0) {
				setG(z, s, 1.0);
			}
			if (s == 1 && k > 0) {
				tmp = exp(-1* u_total [z] / (R * temp));
				//cout << "tmp: " <<tmp<<" ";
				setG(z, s, tmp);
			}
			if (s > 1) {
				setG(z, s, getG_avg(z, s - 1) * getG(z, 1));
			}
			//			if (s < 5) {
			//			cout << "u("<<z<<"):"<< u_total[z]<<" U/(R*temp)="<<u_total[z] / (R * temp);
			//							cout <<" s: "<<s<<" G_avg " <<getG_avg(z,s-1)<<"G(z,s): "<<getG(z,s)<<endl;
			//			}


			//boundary condition
			if (z == zmax - 1) {
				//cout << "JUHU"<<endl;
				setG(z, s, 1.0);
			}

//			if (z==1) {
//				cout << "Mol. :Iteration: "<<k<<" G[" << z << "," << s << "] "<< getG(z,s)<<" u_total["<<z<<"] "<<u_total[z]<<endl;
//			}
		}
	}
}
//G_N
void Molecule::calcG_N(int k, double temp) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = N; s >= 1; s--) {
			int pos = 0;
			pos = s + z * N - N;
			//cout <<"s: "<< s<<endl;
			//initial conditions
			if (s == N && k == 0) {
				setG_N(z, s, 1.0);
			}
			if (s == N && k > 0) {
				tmp = exp(-1* u_total [z] / (R * temp));
				//cout << "tmp: " <<tmp<<" ";
				setG_N(z, s, tmp);
			}
			if (s < N) {
				setG_N(z, s, getG_N_avg(z, s + 1) * getG(z, 1.0));
			}
			//			if (s < 5) {
			//				cout << "u("<<z<<"):"<< u_total[z]<<" U/(R*temp)="<<u_total[z] / (R * temp);
			//								cout <<" s: "<<s<<" G_avg " <<getG_avg(z,s-1)<<"G(z,s): "<<getG(z,s)<<endl;
			//			}


			//boundary condition
			if (z == zmax - 1) {
				//cout << "JUHU"<<endl;
				setG_N(z, s, 1.0);
			}
			//cout << " G[" << z << "," << s << "]=G[" << pos << "] "<< getG(z,s);
			//if ((z*s + 1) % zmax == 0) {

			//}
		}
	}
}

//free chains
void Molecule::calcG_free(int k, double temp) {
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = 1; s <= N; s++) {
			//int pos = 0;
			//pos = s + z * N - N;
			//initial condition
			if (s == 1 && k == 0) {
				setG_free(z, s, 1.0);
			}
			if (s == 1 && k > 0) {
				//cout << "u(z):"<< -u_total[z]<<" -U/(R*temp)="<<-u_total[z] / (R * temp)<<endl;
				//cout << "u(z):"<< -u_total[z]<<endl;
				setG_free(z, s, exp(-u_total[z] / (R * temp)));
			}
			if (s > 1) {
				//cout <<"G_avrg_free " <<getG_free_avg(z,s-1)<<endl;
				setG_free(z, s, getG_free_avg(z, s - 1) * getG(z, 1));
			}
			//boundary condition
			if (z == 1) {
				setG_free(z, s, 0.0);
			}
			//}
		}
	}
}
//calculate fi
void Molecule::calcfi() {
	double tmp;
	total = 0;
	adsorbed_amount = 0;
	theta=0;
	for (int z = 1; z < zmax; z++) {
		fi_z[z] = 0;
		fi_zfree[z] = 0;
		fi_zads[z] = 0;
		for (int s = 1; s <= N; s++) {
			//			tmp = fi_b / ((double) N) * getG(z, s) * getG(z, N - s + 1) / getG(z,
			//					1);
			tmp = fi_b / ((double) N) * getG(z, s) * getG_N(z, s) / getG(z, 1);
			//tmp = fi_b/(double) N;
			//cout << "G(z"<<z<<",s"<<s<<"|1):"<< getG(z, s) << " G(z"<<z<<",s"<<s<<"|N):"<< getG_N(z,s)<<endl;


			//			if (s <5) {
			//			cout <<"z: "<<z<<"fi:" <<tmp;
			//			cout <<" s: "<<s<<" G_avg " <<getG_avg(z,s-1)<<"G(z,s): "<<getG(z,s)<<endl;
			//			}
			setfi(z, s, tmp);
			fi_z[z] = fi_z[z] + tmp;
			//cout << fi_z[z]<<endl;
			//FREE CHAINS
			if (z > 1) {
				tmp = fi_b / (double) N * getG_free(z, s) * getG_free(z, N - s
						+ 1) / getG(z, 1);
				setfi_free(z, s, tmp);

			}

			if (z == 1) {
				tmp = 0.0;
				setfi_free(z, s, 0.0);
			}
			fi_zfree[z] = fi_zfree[z] + tmp;
			fi_zads[z] = fi_z[z] - fi_zfree[z];
		}
		//cout <<" fi_z[z]:"<< fi_z[z]<<endl;
		total = total + fi_z[z];
		adsorbed_amount = adsorbed_amount + fi_zads[z];

	}
}

//G_avg
double Molecule::getG_avg(int z, int s) {
	double result = 0;
	//int pos =0;

	//pos = s + z * N - N;
	//cout << pos;
	//cout << pos +N;
	if (z == 0) {
		exit(1);
	}
	if (z == 1) {
		result = cubic2 * getG(z, s) + cubic1 * getG(z + 1, s);
	}
	if (z > 1 && z < zmax - 1) {
		result = cubic1 * getG(z - 1, s) + cubic2 * getG(z, s) + cubic1 * getG(
				z + 1, s);
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << G[pos+N];
	}
	if (z == zmax - 1) {
		result = cubic1 * getG(z - 1, s) + cubic2 * getG(z, s) + cubic1;
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << cubic1;
	}
	//
	if (z >= zmax) {
		cout << "\nWrong argument for getG_avg(z,s)!";
		exit(1);
	}
	return result;
}

//G_N_avg
double Molecule::getG_N_avg(int z, int s) {
	double result = 0;
	//int pos =0;

	//pos = s + z * N - N;
	//cout << pos;
	//cout << pos +N;
	if (z == 0) {
		exit(1);
	}
	if (z == 1) {
		result = cubic2 * getG_N(z, s) + cubic1 * getG_N(z + 1, s);
	}
	if (z > 1 && z < zmax - 1) {
		result = cubic1 * getG_N(z - 1, s) + cubic2 * getG_N(z, s) + cubic1
				* getG_N(z + 1, s);
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << G[pos+N];
	}
	if (z == zmax - 1) {
		result = cubic1 * getG_N(z - 1, s) + cubic2 * getG_N(z, s) + cubic1;
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << cubic1;
	}
	//
	if (z >= zmax) {
		cout << "\nWrong argument for getG_avg(z,s)!";
		exit(1);
	}
	return result;
}

//G_free_avg
double Molecule::getG_free_avg(int z, int s) {
	double result = 0;
	//int pos =0;
	//pos = s + z * N - N;
	//cout << pos;
	//cout << pos +N;
	if (z == 0) {
		exit(1);
	}
	if (z == 1) {
		result = cubic1 * getG_free(z + 1, s);
	}

	if (z > 1 && z < zmax - 1) {
		result = cubic1 * getG_free(z - 1, s) + cubic2 * getG_free(z, s)
				+ cubic1 * getG_free(z + 1, s);
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << G[pos+N];
	}
	if (z == zmax - 1) {
		result = cubic1 * getG_free(z - 1, s) + cubic2 * getG_free(z, s)
				+ cubic1;
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << cubic1;
	}
	//
	if (z >= zmax) {
		cout << "\nWrong argument for getG_free_avg(z,s)!";
		exit(1);
	}
	return result;
}

//get G (z,s)
//transforms position (z,s) --> pos
double Molecule::getG(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G[pos];
	return result;
}
//get G_N(z,s)
double Molecule::getG_N(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_N[pos];
	return result;
}
//free chains
double Molecule::getG_free(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_free[pos];
	return result;
}

//gets fi(z,s)
double Molecule::getfi(const int z, const int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getfi out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = fi[pos];
	return result;
}

void Molecule::printfiz() {
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		cout << "fi_z[" << z << "]: " << fi_z[z] << endl;
	}
}

void Molecule::setG(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G[pos] = value;
}

void Molecule::setG_N(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_N[pos] = value;
}

void Molecule::setG_free(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_free[pos] = value;
}

void Molecule::setfi(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setfi out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	fi[pos] = value;
}

void Molecule::setfi_free(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setfi_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	fi_free[pos] = value;
}

//fi_avg
double Molecule::getfi_avg(int z) {
	double result = 0;
	if (z == 0) {
		exit(1);
	}
	if (z == 1) {
		result = cubic2 * fi_z[z] + cubic1 * fi_z[z + 1];
	}
	if (z > 1 && z < zmax - 1) {
		result = cubic1 * fi_z[z - 1] + cubic2 * fi_z[z] + cubic1 * fi_z[z + 1];
	}
	if (z == zmax - 1) {
		result = cubic1 * fi_z[z - 1] + cubic2 * fi_z[z] + cubic1 * fi_b;
	}
	return result;
}
//resetting fi_bulk value and related values
void Molecule::setfi_bulk(double fi_copo) {
	fi_b = fi_copo;
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = fi_b;
		fi_zfree[z] = fi_b;
		fi_zads[z] = 0;
		if (z == 0) {
			fi_z[z] = 0.0;
			fi_zfree[z] = 0.0;
		}
	}
	for (int i = 0; i < zmax * N; i++) {
		fi[i] = fi_b / (double) N;
		fi_free[i] = fi_b / (double) N;
	}
	//
}

//add fi_bulk value and related values
void Molecule::addfi_bulk(const double fi_copo) {
	fi_b = fi_b+fi_copo;
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = fi_z[z]+fi_copo;
		fi_zfree[z] = fi_zfree[z]+fi_copo;
		fi_zads[z] = 0;
		if (z == 0) {
			fi_z[z] = 0.0;
			fi_zfree[z] = 0.0;
		}
	}
	for (int i = 0; i < zmax * N; i++) {
		fi[i] = fi_b / (double) N;
		fi_free[i] = fi_b / (double) N;
	}
	//
}

//add fi_z value
void Molecule::addfi_z(const double fi[]) {
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = fi_z[z]+fi[z];
	}
	//
}

//add fi_z value
void Molecule::clearfi_z() {
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = 0;
	}
	//
}

Molecule::~Molecule() {
	// TODO Auto-generated destructor stub
}

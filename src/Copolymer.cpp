/*
 * Copolymer1.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: Chris
 *      //TO DO: G_free is NOT Symmetric: define getG_free_N
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include "Molecule.h"
#include "ScfCalc.h"
#include "Copolymer.h"

Copolymer::Copolymer(int lident, double lfi_b, const string& lPolseq, const Molecule lA,
		const Molecule lB) {
	//Identification
	ident = lident;
	fi_bA = 0;
	fi_bB = 0;
	fi_b = lfi_b;
	N=0;
	cout << "Copolymer#" << ident << " sequence: ";
	//Determining fi_A & fi_B
	//first entry S=0 is empty!!!!!
	//cout<<"LENGTH:"<<lPolseq.length();
	Polseq = new bool[(int) lPolseq.length() + 1];
	for (int s = 0; s < (int) lPolseq.length(); s++) {
		//Polseq[s] = lPolseq[s];
		if (lPolseq[s] == 'A') {
			Polseq[s + 1] = 0;
			cout << "A";
			fi_bA = fi_bA + fi_b / (double) lPolseq.length();
		} else {
			Polseq[s + 1] = 1;
			cout << "B";
			fi_bB = fi_bB + fi_b / (double) lPolseq.length();
		}
		//cout<<N;
		N++;
	}
	//
	for (int z = 0; z < zmax; z++) {
		fi_z[z] = fi_b;
		fi_A[z] = fi_bA;
		fi_B[z] = fi_bB;
		//fi_z[z] = 0.5;
		fi_zfree[z] = fi_b;
		fi_zads[z] = 0;
		fi_zads2[z] = 0;
		fi_tail[z] = 0;
		fi_loop[z] = 0;
		if (z == 0) {
			fi_z[z] = 0.0;
			fi_zfree[z] = 0.0;
		}
		//uA[z]=0;
		//uB[z]=0;
		//alpha[z]=-u_mixing[z];
	}

	cout << endl;
	cout << "Copolymer#" << ident << " length: " << N << endl;
	cout << setprecision(2) << "bulk concentration: " << fi_b << endl
			<< "segment A: " << fi_bA << endl << "segment B: " << fi_bB << endl;
	//setting G;
	//initialisation and declaration of "matrices"
	fi = new double[zmax * N];
	fi_free = new double[zmax * N];
	G = new double[zmax * N];
	G_N = new double[zmax * N];
	G_free = new double[zmax * N];
	G_free_N = new double[zmax * N];
	G_ads = new double[zmax * N];
	G_ads_N = new double[zmax * N];
	for (int i = 0; i < zmax * N; i++) {
		//first entry z=0, s=1 and the last N-1 entries are never used
		fi[i] = fi_b / (double) N;
		fi_free[i] = fi_b / (double) N;
		G[i] = 0.0;
		G_N[i] = 0.0;
		G_free[i] = 0.0;
		G_free_N[i] = 0.0;
		G_ads[i] = 0.0;
		G_ads_N[i] = 0.0;
	}

}
//setting G 1....s
void Copolymer::calcG(int k, double temp, const double uA[], const double uB[]) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = 1; s <= N; s++) {
			//initial conditions
			if (s == 1 && k == 0) {
				setG(z, s, 1.0);
			}
			//next iteration
			//is there a new G(z,s) needed, potential for each segment?No!
			if (Polseq[s] == 0) {
				//cout << "A ";
				tmp = exp(-1* uA [z] / (R * temp));
			} else {
				//cout << "B ";
				tmp = exp(-1* uB [z] / (R * temp));
			}
			if (s == 1 && k > 0) {
				setG(z, s, tmp);
			}
			if (s > 1) {
				setG(z, s, getG_avg(z, s - 1) * tmp);
			}
			//boundary condition
			if (z == zmax - 1) {
				setG(z, s, 1.0);
			}
			//			if (z == 1) {
			//				if (Polseq[s] == 0) {
			//				cout << "Copo. A:Iteration: " << k << " G[" << z << "," << s
			//						<< "] " << getG(z, s) <<" G_avg: " << getG_avg(z, s - 1)<< " u_total[" << z << "] "
			//						<< uA[z] << endl;
			//				} else {
			//					cout << "Copo. B:Iteration: " << k << " G[" << z << "," << s
			//											<< "] " << getG(z, s) <<" G_avg: " << getG_avg(z, s - 1)<< " u_total[" << z << "] "
			//											<< uB[z] << endl;
			//				}
			//			}
		}
	}
}

//setting G N...s
void Copolymer::calcG_N(int k, double temp, const double uA[],
		const double uB[]) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = N; s >= 1; s--) {
			//initial conditions
			if (s == N && k == 0) {
				setG_N(z, s, 1.0);
			}
			//next iteration
			//is there a new G(z,s) needed, potential for each segment?No!
			if (Polseq[s] == 0) {
				//cout << "A ";
				tmp = exp(-1* uA [z] / (R * temp));
			} else {
				//cout << "B ";
				tmp = exp(-1* uB [z] / (R * temp));
			}
			if (s == N && k > 0) {
				setG_N(z, s, tmp);
			}
			if (s < N) {
				setG_N(z, s, getG_N_avg(z, s + 1) * tmp);
			}
			//boundary condition
			if (z == zmax - 1) {
				//cout << "JUHU"<<endl;
				setG_N(z, s, 1.0);
			}
			//			if (z == 1) {
			//							if (Polseq[s] == 0) {
			//							cout << "Copo_N. :Iteration: " << k << " G[" << z << "," << s
			//									<< "] " << getG(z, s) << " u_total[" << z << "] "
			//									<< uA[z] << endl;
			//							} else {
			//								cout << "Copo_N. :Iteration: " << k << " G[" << z << "," << s
			//														<< "] " << getG(z, s) << " u_total[" << z << "] "
			//														<< uB[z] << endl;
			//							}
			//						}

		}
	}
}

//free chains
void Copolymer::calcG_free(int k, double temp, const double uA[],
		const double uB[]) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = 1; s <= N; s++) {
			//int pos = 0;
			//pos = s + z * N - N;
			//initial condition
			if (s == 1 && k == 0) {
				setG_free(z, s, 1.0);
				setG_ads(z, s, 0.0);
			}
			//next iteration
			if (Polseq[s] == 0) {
				//cout << "A ";
				tmp = exp(-1* uA [z] / (R * temp));
			} else {
				//cout << "B ";
				tmp = exp(-1* uB [z] / (R * temp));
			}
			//tmp=getG(z,s);
			if (s == 1 && k > 0) {
				setG_free(z, s, tmp);
			}
			if (s > 1) {
				//cout <<"G_avrg_free " <<getG_free_avg(z,s-1)<<endl;
				setG_free(z, s, getG_free_avg(z, s - 1) * tmp);
			}
			//boundary condition
			if (z == 1) {
				setG_free(z, s, 0.0);
			}
			//setting G_ads = G - G_free
			tmp = getG(z, s) - getG_free(z, s);
			setG_ads(z, s, tmp);
			//cout << "G: "<< getG(z,s)<< " G_free: "<< getG_free(z,s) << " G_ads: "<< getG_ads(z,s)<<endl;

			//}
		}
	}
}

//free chains, no inversion symmetry
void Copolymer::calcG_free_N(int k, double temp, const double uA[],
		const double uB[]) {
	double tmp = 0;
	for (int z = 1; z < zmax; z++) {
		//cout << endl;
		for (int s = N; s >= 1; s--) {
			//int pos = 0;
			//pos = s + z * N - N;
			//initial condition
			if (s == N && k == 0) {
				setG_free_N(z, s, 1.0);
			}
			//next iteration
			if (Polseq[s] == 0) {
				//cout << "A ";
				tmp = exp(-1* uA [z] / (R * temp));
			} else {
				//cout << "B ";
				tmp = exp(-1* uB [z] / (R * temp));
			}
			//tmp=getG(z,s);
			if (s == N && k > 0) {
				setG_free_N(z, s, tmp);
			}
			if (s < N) {
				//cout <<"G_avrg_free " <<getG_free_avg(z,s-1)<<endl;
				setG_free_N(z, s, getG_free_N_avg(z, s + 1) * tmp);
			}
			//boundary condition
			if (z == 1) {
				setG_free_N(z, s, 0.0);
			}
			//setting G_ads = G - G_free
			tmp = getG_N(z, s) - getG_free_N(z, s);
			setG_ads_N(z, s, tmp);
			//cout << "G: "<< getG(z,s)<< " G_free: "<< getG_free(z,s) << " G_ads: "<< getG_ads(z,s)<<endl;

			//}
		}
	}
}

void Copolymer::setG(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G[pos] = value;
	//cout<<value;
}
//get G (z,s)
//transforms position (z,s) --> pos
void Copolymer::setG_N(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_N[pos] = value;
}

void Copolymer::setG_free(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_free[pos] = value;
}

void Copolymer::setG_free_N(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_free_N[pos] = value;
}

void Copolymer::setG_ads(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_ads[pos] = value;
}

void Copolymer::setG_ads_N(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	G_ads_N[pos] = value;
}

double Copolymer::getG(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG out of bounds z:" << z << " s: " << s << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G[pos];
	return result;
}
//get G_N(z,s)
double Copolymer::getG_N(int z, int s) {
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
double Copolymer::getG_free(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG_free out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_free[pos];
	return result;
}

//free chains
double Copolymer::getG_free_N(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG_free out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_free_N[pos];
	return result;
}

//ads chains
double Copolymer::getG_ads(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG_ads out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_ads[pos];
	return result;
}

//ads chains
double Copolymer::getG_ads_N(int z, int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getG_ads_N out of bounds" << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = G_ads_N[pos];
	return result;
}

//get G_average
double Copolymer::getG_avg(int z, int s) {
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
double Copolymer::getG_N_avg(int z, int s) {
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
double Copolymer::getG_free_avg(int z, int s) {
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

//G_free_N_avg
double Copolymer::getG_free_N_avg(int z, int s) {
	double result = 0;
	//int pos =0;
	//pos = s + z * N - N;
	//cout << pos;
	//cout << pos +N;
	if (z == 0) {
		exit(1);
	}
	if (z == 1) {
		result = cubic1 * getG_free_N(z + 1, s);
	}

	if (z > 1 && z < zmax - 1) {
		result = cubic1 * getG_free_N(z - 1, s) + cubic2 * getG_free_N(z, s)
				+ cubic1 * getG_free_N(z + 1, s);
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << G[pos+N];
	}
	if (z == zmax - 1) {
		result = cubic1 * getG_free_N(z - 1, s) + cubic2 * getG_free_N(z, s)
				+ cubic1;
		//cout <<endl<<"G(z,s)="<<G[pos-N]<< " G(z,s)="<<G[pos] << " G(z+1,s)=" << cubic1;
	}
	//
	if (z >= zmax) {
		cout << "\nWrong argument for getG_free_N_avg(z,s)!";
		exit(1);
	}
	return result;
}

//calculate fi
void Copolymer::calcfi(double temp, const double uA[], const double uB[]) {
	double tmp, tmp2;
	total = 0;
	adsorbed_amount = 0;
	loop_amount = 0;
	tail_amount = 0;
	theta = 0;
	for (int z = 1; z < zmax; z++) {
		fi_z[z] = 0;
		fi_A[z] = 0;
		fi_B[z] = 0;
		fi_zfree[z] = 0;
		fi_zads[z] = 0;
		fi_zads2[z] = 0;
		fi_tail[z] = 0;
		fi_loop[z] = 0;
		for (int s = 1; s <= N; s++) {
			//Determine segment type
			if (Polseq[s] == 0) {
				//cout << "A,  s:"<<s<<endl;
				tmp2 = exp(-1* uA [z] / (R * temp));
				tmp = fi_b / ((double) N) * getG(z, s) * getG_N(z, s) / tmp2;
				//fi_A[z] = fi_A[z] + tmp;
			} else {
				//cout << "B ";
				tmp2 = exp(-1* uB [z] / (R * temp));
				tmp = fi_b / ((double) N) * getG(z, s) * getG_N(z, s) / tmp2;
				//fi_B[z] = fi_B[z] + tmp;
			}
			if (z == zmax - 1) {
				tmp = fi_b / ((double) N);
			}
			if (Polseq[s] == 0) {
				fi_A[z] = fi_A[z] + tmp;
			} else {
				fi_B[z] = fi_B[z] + tmp;
			}
			//Normalisation constant??
			setfi(z, s, tmp);
			//fi_z[z] = fi_A[z] + fi_B[z];
			fi_z[z] = fi_z[z] + tmp;
			//free chains
			//TO DO: G_free is NOT Symmetric: define getG_free_N!
			if (z > 1) {
				tmp = fi_b / ((double) N) * getG_free(z, s) * getG_free_N(z, s)
						/ tmp2;
				setfi_free(z, s, tmp);
			}
			if (z == 1) {
				tmp = 0.0;
				setfi_free(z, s, 0.0);
			}
			fi_zfree[z] = fi_zfree[z] + tmp;
			//adsorbed chains: G_ads*G_ads + G_ads*G_free + G_free*G_ads
			//cout << "G: "<< getG(z,s)<< " G_free: "<< getG_free(z,s) << " G_ads: "<< getG_ads(z,s)<<endl;
			tmp = fi_b / ((double) N) * (getG_ads(z, s) * getG_ads_N(z, s)
					+ getG_free(z, s) * getG_ads_N(z, s) + getG_ads(z, s)
					* getG_free_N(z, s)) / tmp2;
			fi_zads2[z] = fi_zads2[z] + tmp;
			//tails: G_ads*G_free + G_free*G_ads
			tmp = fi_b / ((double) N) * (getG_free(z, s) * getG_ads_N(z, s)
					+ getG_ads(z, s) * getG_free_N(z, s)) / tmp2;
			fi_tail[z] = fi_tail[z] + tmp;
			//loops: G_ads*G_ads
			tmp = fi_b / ((double) N) * (getG_ads(z, s) * getG_ads_N(z, s))
					/ tmp2;
			fi_loop[z] = fi_loop[z] + tmp;
			fi_zads[z] = fi_z[z] - fi_zfree[z];
		}
		//cout <<" fi_z[z]:"<< fi_z[z]<<endl;
		total = total + fi_z[z];
		adsorbed_amount = adsorbed_amount + fi_zads[z];
		tail_amount = tail_amount + fi_tail[z];
		loop_amount = loop_amount + fi_loop[z];
	}
}
void Copolymer::setfi(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setG out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	fi[pos] = value;
}

void Copolymer::setfi_free(int z, int s, double value) {
	if (z >= zmax || s > N) {
		cout << "Warning: setfi_free out of bounds" << endl;
		exit(1);
	}
	int pos = 0;
	pos = s + z * N - N;
	fi_free[pos] = value;
}

//gets fi(z,s)
double Copolymer::getfi(const int z, const int s) {
	if (z >= zmax || s > N) {
		cout << "Warning: getfi out of bounds s:" << s << " N:" << N << endl;
		exit(1);
	}
	double result = 0;
	int pos = 0;
	pos = s + z * N - N;
	result = fi[pos];
	return result;
}

//gets fi(z,s)
double Copolymer::getfi_z(int z) {
	if (z >= zmax || z < 1) {
		cout << "Warning: getfi out of bounds z:" << z << endl;
		exit(1);
	}
	return fi_z[z];
}

//partitioning of fi
void Copolymer::partitionfi() {
	for (int z = 1; z < zmax; z++) {
		fi_A[z] = 0;
		fi_B[z] = 0;
		//cout << endl;
		for (int s = 1; s <= N; s++) {
			//cout <<"s:"<<s<<"N: "<<N<<" ";
			if (Polseq[s] == 0) {
				//cout << "fi_A:"<<getfi(s,z)<<endl;
				fi_A[z] = fi_A[z] + getfi(z, s);
			} else {
				//cout <<"fi_B:"<<getfi(s,z)<<endl;
				fi_B[z] = fi_B[z] + getfi(z, s);
			}
		}
		cout << "fi_A[" << z << "]:" << fi_A[z] << endl;
		cout << "fi_B[" << z << "]:" << fi_B[z] << endl;
	}
}

string const Copolymer::getSequence() const {
	string tmp = "A";
	tmp.resize(N);
	for (int s = 1; s <= N; s++) {
		//cout <<"s:"<<s<<"N: "<<N<<" ";
		if (Polseq[s] == 0) {
			tmp[s - 1] = 'A';
		} else {
			tmp[s - 1] = 'B';
		}
	}
	//cout<<tmp.length()<<tmp;
	return tmp;
}
//Calculates probability of chain between z1 and z2...
//double const Copolymer::calcProb(const int z1, const int z2) {
//	double tmp = 0;
//		for (int z = 1; z < zmax; z++) {
//			//cout << endl;
//			for (int s = 1; s <= N; s++) {
//				//initial conditions
//				if (s == 1) {
//					setG(z1, 1, 1.0);
//				} else {
//					setG(z, s, 0.0);
//				}
//				//next iteration s>1
//				//is there a new G(z,s) needed, potential for each segment?No!
//				if (Polseq[s] == 0) {
//					//cout << "A ";
//					tmp = exp(-1* uA [z] / (R * temp));
//				} else {
//					//cout << "B ";
//					tmp = exp(-1* uB [z] / (R * temp));
//				}
//				if (s > 1) {
//					setG(z, s, getG_avg(z, s - 1) * tmp);
//				}
//				//boundary condition
//				if (s == N) {
//					setG(z2, N, 1.0);
//				}
//			}
//		}
//}

Copolymer::~Copolymer() {
	// TODO Auto-generated destructor stub
}

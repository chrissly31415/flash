/*
 * Molecule.h
 *
 *  Created on: Jan 9, 2011
 *      Author: chrissi
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "flash.h"


//using namespace std;

class Molecule {
public:
	Molecule(int lN=1, double lchi_s=0.0,double lchi_12=0.0, double lchi_O=0.0, double lfi_b=0.5);
	virtual ~Molecule();
	//Degree of polymerization
	int N;
	//Flory-Huggins parameters-> array
	double chi_12;
	double chi_O;
	//Surface adhesion parameter (in R*T)
	double chi_s;
	//Concentrations
	double fi_b;
	double* fi;
	double* fi_free;
	double fi_z[zmax];
	double fi_zfree[zmax];
	double fi_zads[zmax];

	//Boltzmann weights
	double* G;
	double* G_N;
	double* G_free;
	//energy parameters
	double u_total[zmax];
	double u_adsorption[zmax];
	double alpha[zmax];
	double alpha_new[zmax];
	double u_mixing[zmax];

	//analysis parameters
	double free;
	double total;
	double excess;
	double adsorbed_amount;
	double theta;
	//functions
	void calcG(int k, double temp);
	void calcG_N(int k, double temp);
	void calcG_free(int k, double temp);
	void calcfi();
	double getG(int z, int s);
	double getG_N(int z, int s);
	double getG_free(int z, int s);
	double getfi(const int z, const int s);
	void printfiz();
	void setG(int z, int s, double value);
	void setG_N(int z, int s, double value);
	void setG_free(int z, int s, double value);
	void setfi(int z, int s, double value);
	void setfi_free(int z, int s, double value);
	double getG_avg(int z, int s);
	double getG_N_avg(int z, int s);
	double getG_free_avg(int z, int s);
	double getfi_avg(int z);
	void setfi_bulk(double fi);
	void addfi_bulk(const double fi);
	void addfi_z(const double fi[]);
	void clearfi_z();
};

#endif /* MOLECULE_H_ */

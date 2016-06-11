/*
 * Copolymer.h
 *
 *  Created on: Jan 31, 2011
 *      Author: jars
 */

#ifndef COPOLYMER_H_
#define COPOLYMER_H_

#include <iostream>
#include <iomanip>
#include "flash.h"
#include "Molecule.h"

using namespace std;

class Copolymer {
public:
	Copolymer(int ident, double lfi_bX, const string& lPolseq,  const Molecule lA,  const Molecule lB);
	virtual ~Copolymer();
	//identification
		int ident;
	//Degree of polymerization, from sequence
		int N;
		//Sequence
		//PolSeqX
		bool* Polseq;
		//define functions print Moleculetype at pos
		//rename Molecule
		Molecule A;
		Molecule B;

		//Concentrations
		double fi_b;
		double fi_bA;
		double fi_bB;
		double* fi;
		double* fi_free;
		double fi_z[zmax];
		double fi_A[zmax];
		double fi_B[zmax];

		double fi_zfree[zmax];
		double fi_zads[zmax];
		double fi_zads2[zmax];
		double fi_tail[zmax];
		double fi_loop[zmax];
		//Boltzmann weights
		double* G;
		double* G_N;
		double* G_free;
		double* G_free_N;
		double* G_ads;
		double* G_ads_N;

		//analysis parameters
		double free;
		double total;
		double excess;
		double theta;
		double adsorbed_amount;
		double loop_amount;
		double tail_amount;
		//functions
		void calcG(int k, double temp,  const double uA[],  const double uB[]);
		void calcG_N(int k, double temp,  const double uA[],  const double uB[]);
		void calcG_free(int k, double temp, const double uA[],  const double uB[]);
		void calcG_free_N(int k, double temp, const double uA[],  const double uB[]);
		void setG(int z, int s, double value);
		void setG_N(int z, int s, double value);
		void setG_free(int z, int s, double value);
		void setG_free_N(int z, int s, double value);
		void setG_ads(int z, int s, double value);
		void setG_ads_N(int z, int s, double value);
		double getG(int z, int s);
		double getG_N(int z, int s);
		double getG_free(int z, int s);
		double getG_free_N(int z, int s);
		double getG_ads(int z, int s);
		double getG_ads_N(int z, int s);
		double getG_avg(int z, int s);
		double getG_N_avg(int z, int s);
		double getG_free_avg(int z, int s);
		double getG_free_N_avg(int z, int s);
		void calcfi(double temp, const double uA[],  const double uB[]);
		void partitionfi();
		double getfi(const int z, const int s);
		void setfi(int z, int s, double value);
		void setfi_free(int z, int s, double value);
		double getfi_z(int z);
		string const getSequence() const;
		//double const calcProb(const int z1, const int z2);
};

#endif /* COPOLYMER_H_ */

//============================================================================
// Name        : flash.cpp
// Author      : christoph loschen
// Version     :
// Copyright   : GNU
// Description : Scheutjens Fleer Theory for Polymer Adsorption
//============================================================================
//TO DO: transfer fi_A, fi_B function, genrerate random sequence, genetic algorithm, correlation length, reading in densities, free energy, radius of gyration, single walk, y-dimension, amino acids, generate random copolymer
#include "flash.h"
#include <iostream>
#include "ScfCalc.h"
#include <math.h>
#include <sys/time.h>

using namespace std;

int main() {
	cout << "FLASH, (c)Christoph Loschen, 2011\n";
	cout
			<< "Lattice based SCF calculation of polymer adsorption at interfaces\n";
	cout
			<< "see e.g. Fleer et al.: Polymer at Interfaces, Chapman & Hall, New York (1993)\n\n";
	cout << "Lattice parameters:\nmax. z-coordinate: " << zmax - 1
			<< "\nmonomer length [nm]: " << mlength;
	cout << "\ncell length [nm]: " << (zmax - 1) * mlength << "\n\n";

	//Prepare timing
	timeval t1, t2;
	gettimeofday(&t1, NULL);
	ScfCalc *myC = new ScfCalc();

	//cout << "\nPARAMETERS:\nLength of polymer 1: " <<(*myC).N_1<<"\n";
	//	cout << "Length of polymer 2: " <<(*myC).N_2<<"\n";
	//	cout << "Bulk concentration polymer 1: "<<(*myC).fi_1b<<"\n";
	//	cout << "Bulk concentration polymer 2: "<<(*myC).fi_2b<<"\n";
	//	cout << "zmax: "<<ScfCalc::zmax<< endl;
	//	(*myC).printfiz((*myC).fi_1z);
	//	(*myC).printfizs((*myC).fi_2,(*myC).N_2);
	//double conca[]={0.2,0.0005,0.001,0.005,0.01,0.02,0.05,0.1};
	//(*myC).startSCFCoPo(40000, 298.15, 0.3);
	for (int i = 0; i < 1; i++) {
		//	conca = conca* i;
		//(*myC).startSCFCoPo(40000, 298.15, conca[i]);
		//(*myC).startSCFGenOpt(40000, 298.15, 0.3);
		//(*myC).startSCFCoPo_Ensemble_old(10000, 298.15, 0.4);
		(*myC).startSCFCoPo_Ensemble(100000, 298.15, 0.2);
		//	(*myC).startSCF(10000, 298.15, conca);
	}
	//Final timing
	gettimeofday(&t2, NULL);
	(*myC).printTiming(t1, t2);
	delete myC;

}

//Constructur
flash::flash() {

}
//Destructor
flash::~flash() {

}


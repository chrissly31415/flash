/*
 * scfCalc.h
 *
 *  Created on: Jan 6, 2011
 *      Author: christoph
 */

#ifndef SCFCALC_H_
#define SCFCALC_H_

#include <iostream>
#include <vector>
#include <random>

#include "Molecule.h"
#include "Copolymer.h"


using namespace std;

class ScfCalc {
public:
	//CONSTRUCTOR
	ScfCalc();
	virtual ~ScfCalc();

	//FUNCTIONS
	//void startSCF(const int iter_max, const double temp, double conca);
	//void startSCF3C(const int iter_max, const double temp, double conca, double concb);
	//void startSCFCoPo(const int iter_max, const double temp, double conca);
	void startSCFGenOpt(const int iter_max, const double temp, double conca);
	void startSCFCoPo_Ensemble(const int iter_max, double temp, double conca);
	//void startSCFCoPo_Ensemble_old(const int iter_max, const double temp, double conca);
	//printing functions, object call with const
	//void printInfo(const Molecule A, const Molecule B, const double temp, bool toFile=false);
	void printInfo3C(const Molecule& A, const Molecule& B, const Molecule& O, const double temp, bool verbose=false, bool toFile=false);
	//void printAdstats(const Molecule A, const Molecule B, const double temp, bool toFile=false);
	//void printAdstats3C(const Molecule A, const Molecule B, const Molecule O, const double temp, bool toFile=false);
	void printAdstatsCopo(const Copolymer X, const Molecule O, const double temp, bool toFile=false);
	void printAdstatsCopo_Ensemble(const std::vector<Copolymer> &lcoposet, const Molecule O, const double all, const double tamount, const double temp, bool toFile=false);
	//void printAdstatsCopo_Ensemble_old(const Copolymer X, const Copolymer Y, const Molecule O, const double temp, bool toFile=false);
	//void printfiz(const Molecule A, const Molecule B, bool toFile=false);
	//void printfiz3C(const Molecule A, const Molecule B, const Molecule O, bool toFile=false);
	void printfizCopo(const Copolymer X, const Molecule O, bool toFile=false);
	void printfizCopo_Ensemble(const std::vector<Copolymer> &lcoposet, const Molecule O, bool toFile=false);
	void printfizCopo_Ensemble3D(const std::vector<Copolymer> &lcoposet);
	//void printfizCopo_Ensemble_old(const Copolymer X, const Copolymer Y, const Molecule O, bool toFile=false);
	void printfis(Copolymer X, bool toStdout=false, bool toFile=true);
	void printTiming(timeval &start, timeval &end);
	void savePotential(const Molecule A, const Molecule B, const Molecule O);
	void setPotential(Copolymer& X,Molecule& A, Molecule& B, Molecule& O, double temp);
	void shuffleSequence(const string oldseq, string& newseq);
	void shuffleBlocks(const string block1, int numb1, const string block2, int numb2, string& newseq);
	void readDist(std::vector<string>& seqset, std::vector<double>& concset);
	void getadsDist(const std::vector<Copolymer>& lcoposet, const std::vector<double>& concset, std::vector<double>& adset);
	void normalizeDist(std::vector<double>& n_set, double nconc);
	void writeDist(const std::vector<Copolymer>& lcoposet, const std::vector<double>& concset);
	void readParams(double& chi_AB, double& chi_As, double& chi_Bs, double& chi_Asol, double& chi_Bsol, double& temp, double& csolvent);

//	void printGzs(double G[], int N);
//	double fi_avg(int z, double fi_z[], double fi_b);
//	double G_avg(int z, int N, double G[]);

	//STATIC VARIABLES
	//static const int zmax=20;
	static const int global_iter_max=1000;
	//static const double temp;


	//VARIABLES
	//int beadtypes;
	double lambda_1[zmax];
	double lambda_2[zmax];
	double lambda_1_start;
	double lambda_2_start;

	//convergence threshhold
	double threshhold;
	double* epsilon;
	double alpha_avg;
	double theta;
	bool converged;
	//old densities
	double fi_zA[zmax];
	double fi_zB[zmax];
	double fi_zO[zmax];
	bool fi_exist;
	double u_A[zmax];
	double u_B[zmax];
	double u_O[zmax];


	//random generator
	//gsl_rng * r;
	std::default_random_engine r;
	//std::random_device rd;
	//std::mt19937 gen(rd());
	std::uniform_real_distribution<double> distribution;

};

#endif /* SCFCALC_H_ */

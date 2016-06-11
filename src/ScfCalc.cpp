/*
 * scfCalc.cpp
 *
 *  Created on: Jan 6, 2011
 *      Author: christoph
 */

#include <boost/regex.hpp>
#include "ScfCalc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>


using namespace std;

//Constructor
ScfCalc::ScfCalc() {
	lambda_1_start = 0.01;
	lambda_2_start = 0.001;
	for (int z = 1; z < zmax; z++) {
		lambda_1[z] = lambda_1_start;
		lambda_2[z] = lambda_2_start;
		//special convergence parameters for surface area to enforce fi=1.0
		if (z < 5) {
			lambda_2[z] = lambda_1_start;
		}

	}
	//	//convergence threshhold
	threshhold = 0.0001;
	fi_exist = false;
	//Random number generator

	//gsl_rng_env_setup();
	//gsl_rng_default_seed = 12663;
	//r = gsl_rng_alloc(gsl_rng_default);
	distribution = std::uniform_real_distribution<double> (0.0, 1.0);

	//gsl_ran_flat(r,0,1);
}
//Standard SCF for single Copolymer
//void ScfCalc::startSCFCoPo(const int iter_max, const double temp, double conca) {
//	double chi_cp_solvent = 0.5;
//	double chi_AB = 0.0;
//	//maximum 3 species : monomer A monomer B, Solvent O
//	//creating polymers N, chi_s, chi_12, chi_O, fi_b
//	Molecule *A = new Molecule(1, 10.0, chi_AB, chi_cp_solvent, 0.0);
//	Molecule *B = new Molecule(1, 0.0, chi_AB, chi_cp_solvent, 0.0);
//	//creating solvent N, chi_s, chi_12, chi_13, fi_b
//	Molecule *O = new Molecule(1, 0.0, chi_cp_solvent, chi_cp_solvent, 0.0);
//	//Copolymer is a new classe: takes Molecule and Sequence
//	string seq = "AAABBBBBBBBBBBABBBBBBBBB";
//	//cout << newseq[1]<< newseq.length();
//	//#Ident, N, fi_b, sequence, A, B
//	Copolymer *X = new Copolymer(1, conca, seq, *A, *B);
//	//resetting correct segment concentrations according to sequence
//	(*A).setfi_bulk((*X).fi_bA);
//	(*B).setfi_bulk((*X).fi_bB);
//	(*O).setfi_bulk(1.0 - (*X).fi_b);
//	//checking old potentials
//	if (fi_exist == true) {
//		cout << "READING OLD POTENTIALS..." << endl;
//		//we have to transfer G(z,s)!
//		//z+1 value has to be set by G_avg!
//		//setPotential(*X, *A, *B, *O, temp);
//	}
//	//information of Molecule A, B, O , temp  verbose, toFile
//	//printInfo3C(*A, *B, *O, temp, false, false);
//
//	//Preparing k loop
//	epsilon = new double[iter_max];
//	double diff = 0.0;
//	//double tmp2=0;
//	//double tmp3=0;
//	converged = false;
//	for (int k = 0; k < iter_max; k++) {
//		epsilon[k] = 0.0;
//		//first iteration && initial guess
//		if (k == 0 && fi_exist == true) {
//			continue;
//		}
//		if (k == 0) {
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z] = (*A).chi_12
//						* ((*B).getfi_avg(z) - (*B).fi_b) * R * temp
//						+ (*A).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				(*B).u_mixing[z] = (*B).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*B).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				(*O).u_mixing[z] = (*O).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*O).chi_O * ((*B).getfi_avg(z) - (*B).fi_b) * R
//								* temp;
//				//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
//				//				u_adsorption[ z]* R * temp;
//				//				(*B).alpha[z] = -1*(*B).u_mixing[z]-(*B).u_adsorption[z]* R * temp;
//				//				(*O).alpha[z] = -1*(*O).u_mixing[z]-(*O).u_adsorption[z]* R * temp;
//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] + (* A ). u_total [z
//				]- (*A).u_adsorption[z]* R * temp;
//				(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//				(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
//				//(*A).alpha[z]=0.1;
//				//(*B).alpha[z]=0.0;
//				}
//			}
//			//update of alpha after first iteration
//			if (k > 0) {
//				for (int z = 1; z < zmax; z++) {
//					(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//					(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z]* R * temp;
//					(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//					(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
//					//cout <<" alpha B:"<<(*B).alpha[z];
//					//New guess for alpha
//					//USE OF AVERAGE ALPHA VALUES!
//					alpha_avg=((*A).alpha[z]+(*B).alpha[z]+(*O).alpha[z])/3;
//					(*A).alpha_new[z]=(*A).alpha[z]+lambda_1[z]*(alpha_avg-(*A).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					(*B).alpha_new[z]=(*B).alpha[z]+lambda_1[z]*(alpha_avg-(*B).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					(*O).alpha_new[z]=(*O).alpha[z]+lambda_1[z]*(alpha_avg-(*O).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					//Error vector
//					diff = sqrt(pow((*A).alpha_new[z]-(*B).alpha_new[z],2)+pow((*A).alpha_new[z]-(*O).alpha_new[z],2)+pow((*B).alpha_new[z]-(*O).alpha_new[z],2));
//					epsilon[k]=epsilon[k] + diff;
//
//					(*A).alpha[z]=(*A).alpha_new[z];
//					(*B).alpha[z]=(*B).alpha_new[z];
//					(*O).alpha[z]=(*O).alpha_new[z];
//					//				if ((*A).alpha_new[z]<0.0) {
//					//					(*A).alpha[z]=-1*(*A).alpha[z];
//					//				}
//
//				}
//			}
//			//convergence test
//			if (sqrt(epsilon[k])<threshhold && k>1000) {
//				cout << fixed << setprecision(4);
//				cout<<"iteration: "<<k<<". Calculation CONVERGED with epsilon<"<< threshhold<<endl;
//				converged = true;
//				break;
//			} else {
//				if (k%500==0) {
//					cout << fixed << setprecision(4);
//					cout << "iteration: " <<k<<" epsilon: "<< sqrt(epsilon[k])<<endl;
//				}
//			}
//			//calculation of u(z)
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//				//cout << (*A).u_mixing[z];
//				//cout << (*A).getfi_avg(z);
//				(*A).u_total[z]=(*A).alpha[z]+(*A).u_mixing[z]+(*A).u_adsorption[z]* R * temp;
//				(*B).u_total[z]=(*B).alpha[z]+(*B).u_mixing[z]+(*B).u_adsorption[z]* R * temp;
//				(*O).u_total[z]=(*O).alpha[z]+(*O).u_mixing[z]+(*O).u_adsorption[z]* R * temp;
//			}
//			//calculation of G(z,s)
//			(*X).calcG(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_N(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_free(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_free_N(k, temp, (*A).u_total, (*B).u_total);
//			//now calculate fi_A and fi_B from fi(s,z)
//			(*X).calcfi(temp, (*A).u_total, (*B).u_total);
//			//(*X).partitionfi();
//
//			//Transfer fi_A from Copolymer to the monomers
//			for (int z = 1; z < zmax; z++) {
//				//cout << "A.fi_z["<<z<<"] "<<(*A).fi_z[z];
//				//with function!
//				(*A).fi_z[z]=(*X).fi_A[z];
//				(*B).fi_z[z]=(*X).fi_B[z];
//				//cout << " A.fi_z["<<z<<"] "<<(*A).fi_z[z]<<endl;
//			}
//			(*O).calcG(k, temp);
//			(*O).calcG_N(k, temp);
//			(*O).calcG_free(k, temp);
//			(*O).calcfi();
//			//Get theta
//			(*X).theta= (*X).fi_z[1] / ((*X).fi_z[1]
//					+ (*O).fi_z[1]);
//			(*O).theta= (*O).fi_z[1] / ((*X).fi_z[1]
//					+ (*O).fi_z[1]);
//			//(*A).printfiz();
//			//(*B).printfiz();
//
//			//(*B).setG(1,1,98.0);
//		}//end of k loop
//		printInfo3C(*A, *B, *O,temp, true, true);
//		printfizCopo(*X, *O, true);
//		printAdstatsCopo(*X, *O, temp, true);
//		printfis(*X, false, true);
//		if (converged != true) {
//			cout<<"\nWARNING: Calculation did NOT CONVERGE after "<<iter_max<<" cycles!"<<endl;
//			fi_exist=false;
//		} //else {
//		//savePotential(*A, *B, *O);
//		//}
//
//		//cout<<(*A).getfi_avg(1);
//
//		delete A;
//		delete B;
//		delete O;
//		delete X;
//	}

void ScfCalc::startSCFGenOpt(const int iter_max, const double temp,
		double conca) {
	double chi_cp_solvent = 0.5;
	double chi_AB = 0.0;
	//maximum 3 species : monomer A monomer B, Solvent O
	//creating polymers N, chi_s, chi_12, chi_O, fi_b
	Molecule *A = new Molecule(1, 10.0, chi_AB, chi_cp_solvent, 0.0);
	Molecule *B = new Molecule(1, 0.0, chi_AB, chi_cp_solvent, 0.0);
	//creating solvent N, chi_s, chi_12, chi_13, fi_b
	Molecule *O = new Molecule(1, 0.0, chi_cp_solvent, chi_cp_solvent, 0.0);
	//Copolymer is a new classe: takes Molecule and Sequence
	//How many sequences are possible? N=24, n1=20, n2=20 -> Nr# = 24!/(20!)(4!)=10626
	//Probability of this sequence: p = pa^4*pB^20=0.5^24=0.00000006
	//(If it were generated by 24 50:50 events!)
	//Probability for a specific sequence * number of specific sequences gives total probability irrespective of order!
	string seq = "AAABBBBBBBBBBBABBBBBBBBB";
	string newseq(seq);
	//shuffleSequence(seq, newseq);
	shuffleBlocks((string) "AAA", 2, (string) "BBBBBBB", 3, newseq);
	//cout << newseq[1]<< newseq.length();
	//#Ident, N, fi_b, sequence, A, B
	Copolymer *X = new Copolymer(1, conca, newseq, *A, *B);
	//resetting correct segment concentrations according to sequence
	(*A).setfi_bulk((*X).fi_bA);
	(*B).setfi_bulk((*X).fi_bB);
	(*O).setfi_bulk(1.0 - (*X).fi_b);
	//checking old potentials
	if (fi_exist == true) {
		cout << "READING OLD POTENTIALS..." << endl;
		//we have to transfer G(z,s)!
		//z+1 value has to be set by G_avg!
		setPotential(*X, *A, *B, *O, temp);
	}
	//information of Molecule A, B, O , temp  verbose, toFile
	//printInfo3C(*A, *B, *O, temp, true, false);

	//Preparing k loop
	epsilon = new double[iter_max];
	double diff = 0.0;
	//double tmp2=0;
	//double tmp3=0;
	converged = false;
	for (int k = 0; k < iter_max; k++) {
		epsilon[k] = 0.0;
		//first iteration && initial guess
		if (k == 0 && fi_exist == true) {
			//continue;
		}
		if (k == 0) {
			for (int z = 1; z < zmax; z++) {
				(*A).u_mixing[z] = (*A).chi_12
						* ((*B).getfi_avg(z) - (*B).fi_b) * R * temp
						+ (*A).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
								* temp;
				//				if (z ==1 ){
				//				cout<<"A.u_mixing: "<<(*A).u_mixing[z] <<"fi_A[z]: "<<(*A).fi_z[z]<<endl;
				//				}
				(*B).u_mixing[z] = (*B).chi_12
						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
						+ (*B).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
								* temp;
				(*O).u_mixing[z] = (*O).chi_12
						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
						+ (*O).chi_O * ((*B).getfi_avg(z) - (*B).fi_b) * R
								* temp;
				//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
				//				u_adsorption[ z]* R * temp;
				//				(*B).alpha[z] = -1*(*B).u_mixing[z]-(*B).u_adsorption[z]* R * temp;
				//				(*O).alpha[z] = -1*(*O).u_mixing[z]-(*O).u_adsorption[z]* R * temp;
				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] + (* A ). u_total [z
				]- (*A).u_adsorption[z]* R * temp;
				(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
				(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
				//(*A).alpha[z]=0.1;
				//(*B).alpha[z]=0.0;
				}
			}
			//update of alpha after first iteration
			if (k > 0) {
				for (int z = 1; z < zmax; z++) {
					(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
					(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
					(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
					(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z]* R * temp;
					(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
					(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
					//cout <<" alpha B:"<<(*B).alpha[z];
					//New guess for alpha
					//USE OF AVERAGE ALPHA VALUES!
					alpha_avg=((*A).alpha[z]+(*B).alpha[z]+(*O).alpha[z])/3.0;
					(*A).alpha_new[z]=(*A).alpha[z]+lambda_1[z]*(alpha_avg-(*A).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
					(*B).alpha_new[z]=(*B).alpha[z]+lambda_1[z]*(alpha_avg-(*B).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
					(*O).alpha_new[z]=(*O).alpha[z]+lambda_1[z]*(alpha_avg-(*O).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
					//Error vector
					diff = sqrt(pow((*A).alpha_new[z]-(*B).alpha_new[z],2)+pow((*A).alpha_new[z]-(*O).alpha_new[z],2)+pow((*B).alpha_new[z]-(*O).alpha_new[z],2));
					epsilon[k]=epsilon[k] + diff;

					(*A).alpha[z]=(*A).alpha_new[z];
					(*B).alpha[z]=(*B).alpha_new[z];
					(*O).alpha[z]=(*O).alpha_new[z];
					//				if ((*A).alpha_new[z]<0.0) {
					//					(*A).alpha[z]=-1*(*A).alpha[z];
					//				}

				}
			}
			//convergence test
			if (sqrt(epsilon[k])<threshhold && k>1000) {
				cout << fixed << setprecision(4);
				cout<<"iteration: "<<k<<". Calculation CONVERGED with epsilon<"<< threshhold<<endl;
				converged = true;
				break;
			} else {

				if (k%500==0) {
					cout << fixed << setprecision(4);
					cout << "iteration: " <<k<<" epsilon: "<< sqrt(epsilon[k])<<endl;
				}
			}
			//calculation of u(z)
			for (int z = 1; z < zmax; z++) {
				(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
				(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
				(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
				//cout << (*A).u_mixing[z];
				//cout << (*A).getfi_avg(z);
				(*A).u_total[z]=(*A).alpha[z]+(*A).u_mixing[z]+(*A).u_adsorption[z]* R * temp;
				(*B).u_total[z]=(*B).alpha[z]+(*B).u_mixing[z]+(*B).u_adsorption[z]* R * temp;
				(*O).u_total[z]=(*O).alpha[z]+(*O).u_mixing[z]+(*O).u_adsorption[z]* R * temp;

			}
			//calculation of G(z,s)
			(*X).calcG(k, temp, (*A).u_total, (*B).u_total);
			(*X).calcG_N(k, temp, (*A).u_total, (*B).u_total);
			(*X).calcG_free(k, temp, (*A).u_total, (*B).u_total);
			(*X).calcG_free_N(k, temp, (*A).u_total, (*B).u_total);
			//now calculate fi_A and fi_B from fi(s,z)
			(*X).calcfi(temp, (*A).u_total, (*B).u_total);
			//cout << (*X).calcProb(2,3);
			//(*X).partitionfi();

			//Transfer fi_A from Copolymer to the monomers
			for (int z = 1; z < zmax; z++) {
				//cout << "A.fi_z["<<z<<"] "<<(*A).fi_z[z];
				//with function!
				(*A).fi_z[z]=(*X).fi_A[z];
				(*B).fi_z[z]=(*X).fi_B[z];
				//cout << " A.fi_z["<<z<<"] "<<(*A).fi_z[z]<<endl;
			}
			(*O).calcG(k, temp);
			(*O).calcG_N(k, temp);
			(*O).calcG_free(k, temp);
			(*O).calcfi();
			//Get theta
			(*X).theta= (*X).fi_z[1] / ((*X).fi_z[1]
					+ (*O).fi_z[1]);
			(*O).theta= (*O).fi_z[1] / ((*X).fi_z[1]
					+ (*O).fi_z[1]);
			//(*A).printfiz();
			//(*B).printfiz();

			//(*B).setG(1,1,98.0);
		}//end of k loop
		//printInfo3C(*A, *B, *O,temp, true, true);
		printfizCopo(*X, *O, false);
		printAdstatsCopo(*X, *O, temp, true);
		printfis(*X, false, true);
		if (converged != true) {
			cout<<"\nWARNING: Calculation did NOT CONVERGE after "<<iter_max<<" cycles!"<<endl;
			fi_exist=false;
		} else {
			savePotential(*A, *B, *O);
		}
		//cout<<(*A).getfi_avg(2);
		//cout<<(*A).getfi_avg(1);
		//		cout<<(*X).getG(zmax-1,1)<<endl;
		//		cout<<(*X).getG(zmax-2,1)<<endl;
		//		cout<<(*X).getG(zmax-3,1)<<endl;
		//		cout<<(*X).getG_avg(zmax-1,1)<<endl;
		//		cout<<(*X).getG_avg(zmax-2,1)<<endl;
		//		cout<<(*X).getG_avg(zmax-3,1)<<endl;
		//		cout<<(*X).getfi_z(zmax-1)<<endl;
		//		cout<<(*O).getG(zmax-1,1)<<endl;
		//		cout<<(*O).getG(zmax-2,1)<<endl;
		delete A;
		delete B;
		delete O;
		delete X;
	}

	//Copo Ensemble, se valgrind
void ScfCalc::startSCFCoPo_Ensemble(const int iter_max, double temp,
		double conca) {
	double tamount = 0;
	double all = 0;
	double chi_Asol = 0.0;
	double chi_Bsol = 0.0;
	double chi_AB = 0.0;
	double chi_As = 0.0;
	double chi_Bs = 0.0;
	double csolvent = 0.0;
	//creating polymers N, chi_s, chi_12, chi_O, fi_b
	readParams(chi_AB, chi_As, chi_Bs, chi_Asol, chi_Bsol, temp, csolvent);
	Molecule *A = new Molecule(1, chi_As, chi_AB, chi_Asol, 0.0);
	Molecule *B = new Molecule(1, chi_Bs, chi_AB, chi_Bsol, 0.0);
	//creating solvent N, chi_s, chi_12, chi_13, fi_b
	Molecule *O = new Molecule(1, 0.0, chi_Asol, chi_Bsol, 0.0);
	//Copolymer is a new classe: takes Molecule and Sequence
	//How many sequences are possible? N=24, n1=20, n2=4 -> Nr# = 24!/(20!)(4!)=10626
	//Probability of this sequence: p = pa^4*pB^20=0.5^24=0.00000006
	//(If it were generated by 24 50:50 events!)
	//Probability for a specific sequence * number of specific sequences gives total probability irrespective of order!
	//-> Binomial Distribution
	std::vector<string> seqset;
	std::vector<double> concset;
	std::vector<double> adset;
	std::vector<double> adset_normal;
	std::vector<string>::iterator itseq;
	std::vector<Copolymer> coposet;
	//reading from file
	int count = 0;
	readDist(seqset, concset);
	normalizeDist(concset, (1.0 - csolvent));
	for (itseq = seqset.begin(); itseq != seqset.end(); ++itseq) {
		//Copolymer *X = new Copolymer(count, concset[count], *itseq, *A, *B);
		//coposet.push_back(*X);
		coposet.push_back(Copolymer(count, concset[count], *itseq, *A, *B));
		count++;
	}

	//Preferential Adsorption of Blockcopolymers at Interfaces as Calculated by Kinetic Monte Carlo Simulations and SCF-Theory
	//Studying INterfacial Weak Boundary Layer Formation with SCF-Theory
	//Vary xi_As versus xi_Bs with length
	//2. Case: Prepolymers i.e. from PU-Reaction, with and without (varying) surfactant
	//	string seq1 = "AAA";
	//	string seq2 = "AAABBBBBBAAA";
	//	string seq3 = "AAABBBBBBBAAABBBBBBBAAA";
	//	string seq4 = "AAABBBBBBBAAABBBBBBBAAABBBBBBBAAABBBBBBBAAA";
	//Findings: Polymer concentrations reduced close to interface, longer chains are depleted
	// ABA dominates in second layer
	// Effect still there if O is  a bit weaker xi_s 10.0 ~9.5 -> chi_s versus theta
	// What happens if O is longer? -> Coverage goes up
	//1. Case: Homopolymers, i.e. polyacrylates and weak boundary formation
	//How to generate a Poisson Distribution
	//a. Count of Strongest adhesive group dominates, then:
	//a. Adhesion promoting group density counts
	//3.Position for adhesion active group
	//string seq1 = "BBBAAABBB";
	//string seq2 = "AAABBBBBB";
	//4.. Case Molecularweight distributions from MC calculation (combination with 1.) for different NCO:OH ratios
	// Effect of monomeric free formulations
	// Effect of 1) solvent concentration, 2) type i.e. bad solvent
	// Effect of endgroups, adhesive or not, comparison with homopolymer etc
	//#Ident, N, fi_b, sequence, A, B
	std::vector<Copolymer>::iterator it;
	for (it = coposet.begin(); it != coposet.end(); ++it) {
		A->addfi_bulk(it->fi_bA);
		B->addfi_bulk(it->fi_bB);
		all = all + it->fi_b;
	}
	cout << "A: c_bulk:" << (*A).fi_b;
	//Handle Monomer/Solvent/Homopolymer separately TO DO: make it flexible
	O->setfi_bulk(1.0 - all);
	//checking old potentials
	//	if (fi_exist == true) {
	//		cout << "READING OLD POTENTIALS..." << endl;
	//		//we have to transfer G(z,s)!
	//		//z+1 value has to be set by G_avg!
	//		setPotential(*X, *A, *B, *O, temp);
	//	}
	//information of Molecule A, B, O , temp  verbose, toFile
	printInfo3C(*A, *B, *O, temp, true, false);
	//Preparing k loop
	epsilon = new double[iter_max];
	double diff = 0.0;
	//double tmp2=0;
	//double tmp3=0;
	converged = false;
	for (int k = 0; k < iter_max; k++) {
		epsilon[k] = 0.0;
		//first iteration && initial guess
		if (k == 0 && fi_exist == true) {
			//continue;
		}
		if (k == 0) {
			for (int z = 1; z < zmax; z++) {
				A->u_mixing[z] = A->chi_12 * (B->getfi_avg(z) - B->fi_b) * R
						* temp + A->chi_O * (O->getfi_avg(z) - O->fi_b) * R
						* temp;
				//				if (z ==1 ){
				//				cout<<"A.u_mixing: "<<A->u_mixing[z] <<"fi_A[z]: "<<A->fi_z[z]<<endl;
				//				}
				B->u_mixing[z] = B->chi_12 * (A->getfi_avg(z) - A->fi_b) * R
						* temp + B->chi_O * (O->getfi_avg(z) - O->fi_b) * R
						* temp;
				O->u_mixing[z] = O->chi_12 * (A->getfi_avg(z) - A->fi_b) * R
						* temp + O->chi_O * (B->getfi_avg(z) - B->fi_b) * R
						* temp;
				//				A->alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
				//				u_adsorption[ z]* R * temp;
				//				B->alpha[z] = -1*B->u_mixing[z]-B->u_adsorption[z]* R * temp;
				//				O->alpha[z] = -1*O->u_mixing[z]-O->u_adsorption[z]* R * temp;
				A->alpha[z] = -1* A ->u_mixing[z] + A->u_total[z]
						- A->u_adsorption[z] * R * temp;
				B->alpha[z] = -1* B ->u_mixing[z] + B->u_total[z]
						- B->u_adsorption[z] * R * temp;
				O->alpha[z] = -1* O ->u_mixing[z] + O->u_total[z]
						- O->u_adsorption[z] * R * temp;
				//A->alpha[z]=0.1;
				//B->alpha[z]=0.0;
			}
		}
		//update of alpha after first iteration
		if (k > 0) {
			for (int z = 1; z < zmax; z++) {
				A->u_mixing[z] = A->chi_12 * (B->getfi_avg(z) - B->fi_b) * R
						* temp + A->chi_O * (O->getfi_avg(z) - O->fi_b) * R
						* temp;
				B->u_mixing[z] = B->chi_12 * (A->getfi_avg(z) - A->fi_b) * R
						* temp + B->chi_O * (O->getfi_avg(z) - O->fi_b) * R
						* temp;
				O->u_mixing[z] = O->chi_12 * (A->getfi_avg(z) - A->fi_b) * R
						* temp + O->chi_O * (B->getfi_avg(z) - B->fi_b) * R
						* temp;
				A->alpha[z] = -1* A ->u_mixing[z] + A->u_total[z]
						- A->u_adsorption[z] * R * temp;
				B->alpha[z] = -1* B ->u_mixing[z] + B->u_total[z]
						- B->u_adsorption[z] * R * temp;
				O->alpha[z] = -1* O ->u_mixing[z] + O->u_total[z]
						- O->u_adsorption[z] * R * temp;
				//cout <<" alpha B:"<<B->alpha[z];
				//New guess for alpha
				//USE OF AVERAGE ALPHA VALUES!
				alpha_avg = (A->alpha[z] + B->alpha[z] + O->alpha[z]) / 3.0;
				A->alpha_new[z] = A->alpha[z] + lambda_1[z] * (alpha_avg
						- A->alpha[z]) - lambda_2[z] * (1.0 - A->fi_z[z]
						- B->fi_z[z] - O->fi_z[z]);
				B->alpha_new[z] = B->alpha[z] + lambda_1[z] * (alpha_avg
						- B->alpha[z]) - lambda_2[z] * (1.0 - A->fi_z[z]
						- B->fi_z[z] - O->fi_z[z]);
				O->alpha_new[z] = O->alpha[z] + lambda_1[z] * (alpha_avg
						- O->alpha[z]) - lambda_2[z] * (1.0 - A->fi_z[z]
						- B->fi_z[z] - O->fi_z[z]);
				//Error vector
				diff = sqrt(pow(A->alpha_new[z] - B->alpha_new[z], 2) + pow(
						A->alpha_new[z] - O->alpha_new[z], 2) + pow(
						B->alpha_new[z] - O->alpha_new[z], 2));
				epsilon[k] = epsilon[k] + diff;

				A->alpha[z] = A->alpha_new[z];
				B->alpha[z] = B->alpha_new[z];
				O->alpha[z] = O->alpha_new[z];
				//				if (A->alpha_new[z]<0.0) {
				//					A->alpha[z]=-1*A->alpha[z];
				//				}

			}
		}
		//convergence test
		if (sqrt(epsilon[k]) < threshhold && k > 1000) {
			cout << fixed << setprecision(4);
			cout << "iteration: " << k
					<< ". Calculation CONVERGED with epsilon<" << threshhold
					<< endl;
			converged = true;
			break;
		} else {

			if (k % 500 == 0) {
				cout << fixed << setprecision(4);
				cout << "iteration: " << k << " epsilon: " << sqrt(epsilon[k])
						<< endl;
			}
		}
		//calculation of u(z)
		for (int z = 1; z < zmax; z++) {
			A->u_mixing[z] = A->chi_12 * (B->getfi_avg(z) - B->fi_b) * R * temp
					+ A->chi_O * (O->getfi_avg(z) - O->fi_b) * R * temp;
			B->u_mixing[z] = B->chi_12 * (A->getfi_avg(z) - A->fi_b) * R * temp
					+ B->chi_O * (O->getfi_avg(z) - O->fi_b) * R * temp;
			O->u_mixing[z] = O->chi_12 * (A->getfi_avg(z) - A->fi_b) * R * temp
					+ O->chi_O * (B->getfi_avg(z) - B->fi_b) * R * temp;
			//cout << A->u_mixing[z];
			//cout << A->getfi_avg(z);
			A->u_total[z] = A->alpha[z] + A->u_mixing[z] + A->u_adsorption[z]
					* R * temp;
			B->u_total[z] = B->alpha[z] + B->u_mixing[z] + B->u_adsorption[z]
					* R * temp;
			O->u_total[z] = O->alpha[z] + O->u_mixing[z] + O->u_adsorption[z]
					* R * temp;

		}
		//clearing fi of segments A&B
		A->clearfi_z();
		B->clearfi_z();
		tamount = 0;
		all = 0;
		//calculation of G(z,s), fi etc for all copolymers
		for (it = coposet.begin(); it != coposet.end(); ++it) {
			//(*it)->fi_b=0.4;
			//resetting correct segment concentrations according to sequence
			it->calcG(k, temp, A->u_total, B->u_total);
			it->calcG_N(k, temp, A->u_total, B->u_total);
			it->calcG_free(k, temp, A->u_total, B->u_total);
			it->calcG_free_N(k, temp, A->u_total, B->u_total);
			it->calcfi(temp, A->u_total, B->u_total);
			//Transfer fi_A from Copolymer to the monomers
			A->addfi_z(it->fi_A);
			B->addfi_z(it->fi_B);
			tamount = tamount + it->fi_z[1];
			all = all + it->total;
			//cout << (*it)->fi_b;
		}
		//Solvent/Monomer separately
		O->calcG(k, temp);
		O->calcG_N(k, temp);
		O->calcG_free(k, temp);
		O->calcfi();
		//get connected properties
		tamount = tamount + O->fi_z[1];
		all = all + O->total;
		//Get theta
		for (it = coposet.begin(); it != coposet.end(); ++it) {
			it->theta = it->fi_z[1] / tamount;
		}
		O->theta = O->fi_z[1] / tamount;
	}//end of k loop
	printInfo3C(*A, *B, *O, temp, true, true);
	printfizCopo_Ensemble(coposet, *O, true);
	printfizCopo_Ensemble3D(coposet);
	printAdstatsCopo_Ensemble(coposet, *O, all, tamount, temp, true);
	getadsDist(coposet, concset, adset);
	writeDist(coposet, adset);
	//normalizeDist(adset);
	//printfis(*X, false, true);
	if (converged != true) {
		cout << "\nWARNING: Calculation did NOT CONVERGE after " << iter_max
				<< " cycles!" << endl;
		fi_exist = false;
	} else {
		savePotential(*A, *B, *O);
	}
	delete A;
	delete B;
	delete O;
	//	delete X1;
	//	delete X2;
	//	delete X3;
	//	delete X4;
	//	delete X5;
	//	delete X6;
	//	delete X7;
	//	delete X8;
	//	delete X9;
	//	delete X10;
	//delete coposet2;
}

//void ScfCalc::startSCFCoPo_Ensemble_old(const int iter_max, const double temp,
//		double conca) {
//	double chi_cp_solvent = 0.0;
//	double chi_AB = 0.0;
//	//maximum 3 species : monomer A monomer B, Solvent O
//	//creating polymers N, chi_s, chi_12, chi_O, fi_b
//	Molecule *A = new Molecule(1, 5.0, chi_AB, chi_cp_solvent, 0.0);
//	Molecule *B = new Molecule(1, 0.0, chi_AB, chi_cp_solvent, 0.0);
//	//creating solvent N, chi_s, chi_12, chi_13, fi_b
//	Molecule *O = new Molecule(1, 0.0, chi_cp_solvent, chi_cp_solvent, 0.0);
//	//Copolymer is a new classe: takes Molecule and Sequence
//	//How many sequences are possible? N=24, n1=20, n2=4 -> Nr# = 24!/(20!)(4!)=10626
//	//Probability of this sequence: p = pa^4*pB^20=0.5^24=0.00000006
//	//(If it were generated by 24 50:50 events!)
//	//Probability for a specific sequence * number of specific sequences gives total probability irrespective of order!
//	//-> Binomial Distribution
//
//	//Preferential Adsorption of Blockcopolymers at Interfaces as Calculated by Kinetic Monte Carlo Simulations and SCF-Theory
//	string seq1 = "AAAAABBBBBAAAAA";
//	string seq2 = "ABBBBB";
//	//string newseq(seq);
//	//shuffleSequence(seq, newseq);
//	//shuffleBlocks((string) "AAA", 2, (string) "BBBBBBB", 3, newseq);
//	//cout << newseq[1]<< newseq.length();
//	//#Ident, N, fi_b, sequence, A, B
//	Copolymer *X = new Copolymer(1, conca, seq1, *A, *B);
//	Copolymer *Y = new Copolymer(2, conca, seq2, *A, *B);
//	//resetting correct segment concentrations according to sequence
//	(*A).setfi_bulk((*X).fi_bA + (*Y).fi_bA);
//	(*B).setfi_bulk((*X).fi_bB + (*Y).fi_bB);
//	(*O).setfi_bulk(1.0 - (*X).fi_b - (*Y).fi_b);
//	//checking old potentials
//	//	if (fi_exist == true) {
//	//		cout << "READING OLD POTENTIALS..." << endl;
//	//		//we have to transfer G(z,s)!
//	//		//z+1 value has to be set by G_avg!
//	//		setPotential(*X, *A, *B, *O, temp);
//	//	}
//	//information of Molecule A, B, O , temp  verbose, toFile
//	printInfo3C(*A, *B, *O, temp, true, false);
//
//	//Preparing k loop
//	epsilon = new double[iter_max];
//	double diff = 0.0;
//	//double tmp2=0;
//	//double tmp3=0;
//	converged = false;
//	for (int k = 0; k < iter_max; k++) {
//		epsilon[k] = 0.0;
//		//first iteration && initial guess
//		if (k == 0 && fi_exist == true) {
//			//continue;
//		}
//		if (k == 0) {
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z] = (*A).chi_12
//						* ((*B).getfi_avg(z) - (*B).fi_b) * R * temp
//						+ (*A).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				//				if (z ==1 ){
//				//				cout<<"A.u_mixing: "<<(*A).u_mixing[z] <<"fi_A[z]: "<<(*A).fi_z[z]<<endl;
//				//				}
//				(*B).u_mixing[z] = (*B).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*B).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				(*O).u_mixing[z] = (*O).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*O).chi_O * ((*B).getfi_avg(z) - (*B).fi_b) * R
//								* temp;
//				//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
//				//				u_adsorption[ z]* R * temp;
//				//				(*B).alpha[z] = -1*(*B).u_mixing[z]-(*B).u_adsorption[z]* R * temp;
//				//				(*O).alpha[z] = -1*(*O).u_mixing[z]-(*O).u_adsorption[z]* R * temp;
//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] + (* A ). u_total [z
//				]- (*A).u_adsorption[z]* R * temp;
//				(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//				(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
//				//(*A).alpha[z]=0.1;
//				//(*B).alpha[z]=0.0;
//				}
//			}
//			//update of alpha after first iteration
//			if (k > 0) {
//				for (int z = 1; z < zmax; z++) {
//					(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//					(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z]* R * temp;
//					(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//					(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
//					//cout <<" alpha B:"<<(*B).alpha[z];
//					//New guess for alpha
//					//USE OF AVERAGE ALPHA VALUES!
//					alpha_avg=((*A).alpha[z]+(*B).alpha[z]+(*O).alpha[z])/3.0;
//					(*A).alpha_new[z]=(*A).alpha[z]+lambda_1[z]*(alpha_avg-(*A).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					(*B).alpha_new[z]=(*B).alpha[z]+lambda_1[z]*(alpha_avg-(*B).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					(*O).alpha_new[z]=(*O).alpha[z]+lambda_1[z]*(alpha_avg-(*O).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					//Error vector
//					diff = sqrt(pow((*A).alpha_new[z]-(*B).alpha_new[z],2)+pow((*A).alpha_new[z]-(*O).alpha_new[z],2)+pow((*B).alpha_new[z]-(*O).alpha_new[z],2));
//					epsilon[k]=epsilon[k] + diff;
//
//					(*A).alpha[z]=(*A).alpha_new[z];
//					(*B).alpha[z]=(*B).alpha_new[z];
//					(*O).alpha[z]=(*O).alpha_new[z];
//					//				if ((*A).alpha_new[z]<0.0) {
//					//					(*A).alpha[z]=-1*(*A).alpha[z];
//					//				}
//
//				}
//			}
//			//convergence test
//			if (sqrt(epsilon[k])<threshhold && k>1000) {
//				cout << fixed << setprecision(4);
//				cout<<"iteration: "<<k<<". Calculation CONVERGED with epsilon<"<< threshhold<<endl;
//				converged = true;
//				break;
//			} else {
//
//				if (k%500==0) {
//					cout << fixed << setprecision(4);
//					cout << "iteration: " <<k<<" epsilon: "<< sqrt(epsilon[k])<<endl;
//				}
//			}
//			//calculation of u(z)
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//				//cout << (*A).u_mixing[z];
//				//cout << (*A).getfi_avg(z);
//				(*A).u_total[z]=(*A).alpha[z]+(*A).u_mixing[z]+(*A).u_adsorption[z]* R * temp;
//				(*B).u_total[z]=(*B).alpha[z]+(*B).u_mixing[z]+(*B).u_adsorption[z]* R * temp;
//				(*O).u_total[z]=(*O).alpha[z]+(*O).u_mixing[z]+(*O).u_adsorption[z]* R * temp;
//
//			}
//			//calculation of G(z,s)
//			(*X).calcG(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_N(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_free(k, temp, (*A).u_total, (*B).u_total);
//			(*X).calcG_free_N(k, temp, (*A).u_total, (*B).u_total);
//			//now calculate fi_A and fi_B from fi(s,z)
//			(*X).calcfi(temp, (*A).u_total, (*B).u_total);
//			//now second polymer
//			(*Y).calcG(k, temp, (*A).u_total, (*B).u_total);
//			(*Y).calcG_N(k, temp, (*A).u_total, (*B).u_total);
//			(*Y).calcG_free(k, temp, (*A).u_total, (*B).u_total);
//			(*Y).calcG_free_N(k, temp, (*A).u_total, (*B).u_total);
//			//now calculate fi_A and fi_B from fi(s,z)
//			(*Y).calcfi(temp, (*A).u_total, (*B).u_total);
//
//			//cout << (*X).calcProb(2,3);
//			//(*X).partitionfi();
//
//			//Transfer fi_A from Copolymer to the monomers
//			for (int z = 1; z < zmax; z++) {
//				//cout << "A.fi_z["<<z<<"] "<<(*A).fi_z[z];
//				//with function!
//				(*A).fi_z[z]=(*X).fi_A[z]+(*Y).fi_A[z];
//				(*B).fi_z[z]=(*X).fi_B[z]+(*Y).fi_B[z];
//				//cout << " A.fi_z["<<z<<"] "<<(*A).fi_z[z]<<endl;
//			}
//			(*O).calcG(k, temp);
//			(*O).calcG_N(k, temp);
//			(*O).calcG_free(k, temp);
//			(*O).calcfi();
//			//Get theta
//			(*X).theta= (*X).fi_z[1] / ((*X).fi_z[1]+(*Y).fi_z[1]
//					+ (*O).fi_z[1]);
//			(*Y).theta= (*Y).fi_z[1] / ((*X).fi_z[1]+(*Y).fi_z[1]
//					+ (*O).fi_z[1]);
//			(*O).theta= (*O).fi_z[1] / ((*X).fi_z[1]+(*Y).fi_z[1]
//					+ (*O).fi_z[1]);
//			//(*A).printfiz();
//			//(*B).printfiz();
//
//			//(*B).setG(1,1,98.0);
//		}//end of k loop
//		printInfo3C(*A, *B, *O,temp, true, true);
//		printfizCopo_Ensemble_old(*X, *Y, *O, true);
//		printAdstatsCopo_Ensemble_old(*X,*Y, *O, temp, true);
//		printfis(*X, false, true);
//		if (converged != true) {
//			cout<<"\nWARNING: Calculation did NOT CONVERGE after "<<iter_max<<" cycles!"<<endl;
//			fi_exist=false;
//		} else {
//			savePotential(*A, *B, *O);
//		}
//		//cout<<(*A).getfi_avg(2);
//		//cout<<(*A).getfi_avg(1);
//		//		cout<<(*X).getG(zmax-1,1)<<endl;
//		//		cout<<(*X).getG(zmax-2,1)<<endl;
//		//		cout<<(*X).getG(zmax-3,1)<<endl;
//		//		cout<<(*X).getG_avg(zmax-1,1)<<endl;
//		//		cout<<(*X).getG_avg(zmax-2,1)<<endl;
//		//		cout<<(*X).getG_avg(zmax-3,1)<<endl;
//		//		cout<<(*X).getfi_z(zmax-1)<<endl;
//		//		cout<<(*O).getG(zmax-1,1)<<endl;
//		//		cout<<(*O).getG(zmax-2,1)<<endl;
//		delete A;
//		delete B;
//		delete O;
//		delete X;
//	}
//
//void ScfCalc::startSCF(const int iter_max, const double temp, double conca) {
//	//creating polymers N, chi_s, chi_12, chi_O, fi_b
//	//maximum 3 species : monomer A monomer B, Solvent O
//	Molecule *A = new Molecule(100, 2.0, 0.25, 0.0, conca);
//	Molecule *B = new Molecule(2, 0.0, 0.25, 0.0, 1.0 - conca);
//	//creating solvent N, chi_s, chi_1, chi_2, fi_b
//	//USE OF AVERAGE VALUES!!
//	//Molecule *O = new Molecule(1, 0.0, 0.0, 0.0, 1.0 - conca);
//	//information of Molecule A, B , print to File
//	printInfo(*A, *B, false);
//	//beginning of k loop
//	epsilon = new double[iter_max];
//	double diff = 0.0;
//	//double tmp2=0;
//	//double tmp3=0;
//	for (int k = 0; k < iter_max; k++) {
//		epsilon[k] = 0.0;
//		//first iteration
//		if (k == 0) {
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z] = (*A).chi_12
//						* ((*B).getfi_avg(z) - (*B).fi_b) * R * temp;
//				(*B).u_mixing[z] = (*B).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp;
//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
//				u_adsorption[ z]* R * temp+0.001;
//				(*B).alpha[z] = -1*(*B).u_mixing[z]-(*B).u_adsorption[z]* R * temp;
//				//(*A).alpha[z]=0.1;
//				//(*B).alpha[z]=0.0;
//				}
//			}
//			//update of alpha after first iteration
//			if (k > 0) {
//				for (int z = 1; z < zmax; z++) {
//					(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//					(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp;
//					(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z]* R * temp;
//					(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//					//(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z];
//					//(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z];
//
//					//cout <<"alpha A:"<<(*A).alpha[z];
//					//cout <<" alpha B:"<<(*B).alpha[z];
//					//New guess for alpha
//					//USE OF AVERAGE ALPHA VALUES!
//					alpha_avg=((*A).alpha[z]+(*B).alpha[z])/2;
//					(*A).alpha_new[z]=(*A).alpha[z]+lambda_1[z]*(alpha_avg-(*A).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]);
//					//(*A).alpha_new[z]=(*A).alpha[z]+0.5*((*B).alpha[z]-(*A).alpha[z]);
//					//cout << (*B).fi_z[z] <<endl;
//					(*B).alpha_new[z]=(*B).alpha[z]+lambda_1[z]*(alpha_avg-(*B).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]);
//					//(*B).alpha_new[z]=(*B).alpha[z]+0.5*lambda_1*((*A).alpha[z]-(*B).alpha[z]);
//					diff = pow((*A).alpha_new[z]-(*B).alpha_new[z],2);
//					//				if (diff > threshhold) {
//					//					cout <<z<< " WARNING: diff: "<<diff<<" > "<< threshhold<<endl;
//					//					lambda_1[z]=lambda_1_start*5.0;
//					//				} else {
//					//					lambda_1[z]=lambda_1_start;
//					//				}
//
//					//if (k == 999) {
//					//cout << z << setprecision(5)<< " "<<diff<< endl;
//					//}
//					epsilon[k]=epsilon[k] + diff;
//					//if (z < 3) {
//
//					//					cout << "u(z)"<< (*A).u_total[z] << " alpha:"<<(*A).alpha[z]<<" u_mixing"<<(*A).u_mixing[z]<<" u_ads:"<<(*A).u_adsorption[z]<<endl;
//					//					cout << "u(z)"<< (*B).u_total[z] << " alpha:"<<(*B).alpha[z]<< " u_mixing"<<(*B).u_mixing[z]<<" u_ads:"<<(*B).u_adsorption[z]<<endl;
//					//					cout <<" alpha_new A:"<<(*A).alpha_new[z]<< " alpha A:" << (*A).alpha[z]<<endl;
//					//					cout <<" alpha_new B:"<<(*B).alpha_new[z]<< " alpha B:" << (*B).alpha[z]<<endl;
//					//}
//					(*A).alpha[z]=(*A).alpha_new[z];
//					(*B).alpha[z]=(*B).alpha_new[z];
//					//				if ((*A).alpha_new[z]<0.0) {
//					//					(*A).alpha[z]=-1*(*A).alpha[z];
//					//				}
//
//				}
//			}
//			//convergence test
//			if (sqrt(epsilon[k])<threshhold && k>1000) {
//
//				cout << fixed << setprecision(4);
//				cout<<"iteration: "<<k<<". Calculation CONVERGED with epsilon<"<< threshhold<<endl;
//
//				break;
//			} else {
//				if (k%500==0) {
//					cout << fixed << setprecision(4);
//					cout << "iteration: " <<k<<" epsilon: "<< sqrt(epsilon[k])<<endl;
//				}
//			}
//			//calculation of u(z)
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//				(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp;
//				//cout << (*A).u_mixing[z];
//				//cout << (*A).getfi_avg(z);
//				(*A).u_total[z]=(*A).alpha[z]+(*A).u_mixing[z]+(*A).u_adsorption[z]* R * temp;
//				(*B).u_total[z]=(*B).alpha[z]+(*B).u_mixing[z]+(*B).u_adsorption[z]* R * temp;
//				//cout << "u(z)"<< (*A).u_total[z] << " alpha:"<<(*A).alpha[z]<<" u_mixing"<<(*A).u_mixing[z]<<" u_ads:"<<(*A).u_adsorption[z]<<endl;
//				//cout << "u(z)"<< (*B).u_total[z] << " alpha:"<<(*B).alpha[z]<< " u_mixing"<<(*B).u_mixing[z]<<" u_ads:"<<(*B).u_adsorption[z]<<endl;
//			}
//			//calculation of G(z,s)
//			(*A).calcG(k, temp);
//			(*B).calcG(k, temp);
//			(*A).calcG_N(k, temp);
//			(*B).calcG_N(k, temp);
//			(*A).calcG_free(k, temp);
//			(*B).calcG_free(k, temp);
//			(*A).calcfi();
//			(*B).calcfi();
//			//(*A).printfiz();
//			//(*B).printfiz();
//
//			//(*B).setG(1,1,98.0);
//			//cout <<endl<< (*B).getG(1,1);
//
//		}//end of k loop
//		printInfo(*A, *B, temp, true);
//		printfiz(*A, *B, true);
//		printAdstats(*A, *B, temp, true);
//		//cout<<(*A).getfi_avg(2);
//		//cout<<(*A).getfi_avg(1);
//		delete A;
//		delete B;
//	}
//
//void ScfCalc::startSCF3C(const int iter_max, const double temp, double conca,
//		double concb) {
//	//creating polymers N, chi_s, chi_12, chi_O, fi_b
//	//maximum 3 species : monomer A monomer B, Solvent O
//	Molecule *A = new Molecule(8, 2.0, 0.0, 0.29, conca);
//	Molecule *B = new Molecule(1, 1.0, 0.0, 0.29, concb);
//	//creating solvent N, chi_s, chi_12, chi_13, fi_b
//	Molecule *O = new Molecule(1, 0.0, 0.29, 0.29, 1.0 - conca - concb);
//	//information of Molecule A, B , print to File
//	printInfo3C(*A, *B, *O, false, false);
//	//beginning of k loop
//	epsilon = new double[iter_max];
//	double diff = 0.0;
//	//double tmp2=0;
//	//double tmp3=0;
//	for (int k = 0; k < iter_max; k++) {
//		epsilon[k] = 0.0;
//		//first iteration
//		if (k == 0) {
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z] = (*A).chi_12
//						* ((*B).getfi_avg(z) - (*B).fi_b) * R * temp
//						+ (*A).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				(*B).u_mixing[z] = (*B).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*B).chi_O * ((*O).getfi_avg(z) - (*O).fi_b) * R
//								* temp;
//				(*O).u_mixing[z] = (*O).chi_12
//						* ((*A).getfi_avg(z) - (*A).fi_b) * R * temp
//						+ (*O).chi_O * ((*B).getfi_avg(z) - (*B).fi_b) * R
//								* temp;
//				;
//				(*A).alpha[z] = -1* (* A ) .u_mixing [ z ] - (* A ).
//				u_adsorption[ z]* R * temp+0.001;
//				(*B).alpha[z] = -1*(*B).u_mixing[z]-(*B).u_adsorption[z]* R * temp;
//				(*O).alpha[z] = -1*(*O).u_mixing[z]-(*O).u_adsorption[z]* R * temp;
//				//(*A).alpha[z]=0.1;
//				//(*B).alpha[z]=0.0;
//				}
//			}
//			//update of alpha after first iteration
//			if (k > 0) {
//				for (int z = 1; z < zmax; z++) {
//					(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//					(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//
//					(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z]* R * temp;
//					(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z]* R * temp;
//					(*O).alpha[z]=-1*(*O).u_mixing[z]+(*O).u_total[z]-(*O).u_adsorption[z]* R * temp;
//					//					if (z==10&& k%500==0) {
//					//						//cout << "A.u_mixing:"<<(*A).u_mixing[z]<< endl;
//					//						//cout << "B.u_mixing:"<<(*B).u_mixing[z]<< endl;
//					//						cout<< "B.getfi_avg: "<< (*B).getfi_avg(z) << " B.fi_B: "<<(*B).fi_b;
//					//						cout <<" "<<(*B).fi_z[z-2]<<" "<<(*B).fi_z[z-1]<<" "<<(*B).fi_z[z]<<endl;
//					//						cout << "O.u_mixing:"<<(*O).u_mixing[z]<< endl;
//					//					}
//
//					//(*B).alpha[z]=-1*(*B).u_mixing[z]+(*B).u_total[z]-(*B).u_adsorption[z];
//					//(*A).alpha[z]=-1*(*A).u_mixing[z]+(*A).u_total[z]-(*A).u_adsorption[z];
//
//					//cout <<"alpha A:"<<(*A).alpha[z];
//					//cout <<" alpha B:"<<(*B).alpha[z];
//					//New guess for alpha
//					//USE OF AVERAGE ALPHA VALUES!
//					alpha_avg=((*A).alpha[z]+(*B).alpha[z]+(*O).alpha[z])/3;
//					(*A).alpha_new[z]=(*A).alpha[z]+lambda_1[z]*(alpha_avg-(*A).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					//(*A).alpha_new[z]=(*A).alpha[z]+0.5*((*B).alpha[z]-(*A).alpha[z]);
//					//cout << (*B).fi_z[z] <<endl;
//					(*B).alpha_new[z]=(*B).alpha[z]+lambda_1[z]*(alpha_avg-(*B).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//					//(*B).alpha_new[z]=(*B).alpha[z]+0.5*lambda_1*((*A).alpha[z]-(*B).alpha[z]);
//					(*O).alpha_new[z]=(*O).alpha[z]+lambda_1[z]*(alpha_avg-(*O).alpha[z])-lambda_2[z]*(1.0-(*A).fi_z[z]-(*B).fi_z[z]-(*O).fi_z[z]);
//
//					diff = sqrt(pow((*A).alpha_new[z]-(*B).alpha_new[z],2)+pow((*A).alpha_new[z]-(*O).alpha_new[z],2)+pow((*B).alpha_new[z]-(*O).alpha_new[z],2));
//					//				if (diff > threshhold) {
//					//					cout <<z<< " WARNING: diff: "<<diff<<" > "<< threshhold<<endl;
//					//					lambda_1[z]=lambda_1_start*5.0;
//					//				} else {
//					//					lambda_1[z]=lambda_1_start;
//					//				}
//
//					//if (k == 999) {
//					//cout << z << setprecision(5)<< " "<<diff<< endl;
//					//}
//					epsilon[k]=epsilon[k] + diff;
//					//if (z < 3) {
//
//					//					cout << "u(z)"<< (*A).u_total[z] << " alpha:"<<(*A).alpha[z]<<" u_mixing"<<(*A).u_mixing[z]<<" u_ads:"<<(*A).u_adsorption[z]<<endl;
//					//					cout << "u(z)"<< (*B).u_total[z] << " alpha:"<<(*B).alpha[z]<< " u_mixing"<<(*B).u_mixing[z]<<" u_ads:"<<(*B).u_adsorption[z]<<endl;
//					//					cout <<" alpha_new A:"<<(*A).alpha_new[z]<< " alpha A:" << (*A).alpha[z]<<endl;
//					//					cout <<" alpha_new B:"<<(*B).alpha_new[z]<< " alpha B:" << (*B).alpha[z]<<endl;
//					//}
//					(*A).alpha[z]=(*A).alpha_new[z];
//					(*B).alpha[z]=(*B).alpha_new[z];
//					(*O).alpha[z]=(*O).alpha_new[z];
//					//				if ((*A).alpha_new[z]<0.0) {
//					//					(*A).alpha[z]=-1*(*A).alpha[z];
//					//				}
//
//				}
//			}
//			//convergence test
//			if (sqrt(epsilon[k])<threshhold && k>1000) {
//
//				cout << fixed << setprecision(4);
//				cout<<"iteration: "<<k<<". Calculation CONVERGED with epsilon<"<< threshhold<<endl;
//
//				break;
//			} else {
//				if (k%500==0) {
//					cout << fixed << setprecision(4);
//					//cout <<endl<<"G(z=1,1)"<< (*B).getG(1,1) <<" G(z=1,N) "<< (*B).getG(1,(*B).N)<<endl;
//					//cout <<"fi(z=1,1)"<< (*B).getfi(1,1) <<" fi(z=1,N-2) "<< (*B).getfi(1,(*B).N-2)<<endl;
//					//	cout <<endl<<"fi_B(z=1)"<< (*B).fi_z[1]<<"fi_B(z=10)"<< (*B).fi_z[10] <<endl;
//					//cout << " alpha: "<<alpha_avg<<endl;
//					cout << "iteration: " <<k<<" epsilon: "<< sqrt(epsilon[k])<<endl;
//				}
//			}
//			//calculation of u(z)
//			for (int z = 1; z < zmax; z++) {
//				(*A).u_mixing[z]=(*A).chi_12*((*B).getfi_avg(z)-(*B).fi_b)*R*temp+(*A).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*B).u_mixing[z]=(*B).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*B).chi_O*((*O).getfi_avg(z)-(*O).fi_b)*R*temp;
//				(*O).u_mixing[z]=(*O).chi_12*((*A).getfi_avg(z)-(*A).fi_b)*R*temp+(*O).chi_O*((*B).getfi_avg(z)-(*B).fi_b)*R*temp;
//				//cout << (*A).u_mixing[z];
//				//cout << (*A).getfi_avg(z);
//				(*A).u_total[z]=(*A).alpha[z]+(*A).u_mixing[z]+(*A).u_adsorption[z]* R * temp;
//				(*B).u_total[z]=(*B).alpha[z]+(*B).u_mixing[z]+(*B).u_adsorption[z]* R * temp;
//				(*O).u_total[z]=(*O).alpha[z]+(*O).u_mixing[z]+(*O).u_adsorption[z]* R * temp;
//				//cout << "u(z)"<< (*A).u_total[z] << " alpha:"<<(*A).alpha[z]<<" u_mixing"<<(*A).u_mixing[z]<<" u_ads:"<<(*A).u_adsorption[z]<<endl;
//				//cout << "u(z)"<< (*B).u_total[z] << " alpha:"<<(*B).alpha[z]<< " u_mixing"<<(*B).u_mixing[z]<<" u_ads:"<<(*B).u_adsorption[z]<<endl;
//			}
//			//calculation of G(z,s)
//			(*A).calcG(k, temp);
//			(*B).calcG(k, temp);
//			(*O).calcG(k, temp);
//			(*A).calcG_N(k, temp);
//			(*B).calcG_N(k, temp);
//			(*O).calcG_N(k, temp);
//			(*A).calcG_free(k, temp);
//			(*B).calcG_free(k, temp);
//			(*O).calcG_free(k, temp);
//			(*A).calcfi();
//			(*B).calcfi();
//			(*O).calcfi();
//			//(*A).printfiz();
//			//(*B).printfiz();
//
//			//(*B).setG(1,1,98.0);
//
//
//		}//end of k loop
//		printInfo3C(*A, *B, *O,temp, true, true);
//		//printfiz3C(*A, *B, *O, true);
//		printAdstats3C(*A, *B, *O, temp, true);
//		//cout<<(*A).getfi_avg(2);
//		//cout<<(*A).getfi_avg(1);
//		delete A;
//		delete B;
//		delete O;
//	}

//void ScfCalc::printInfo(const Molecule A, const Molecule B, const double temp,
//		bool toFile) {
//	cout << fixed << setprecision(2);
//	cout << endl << "Printing general info:" << endl;
//	cout << "Species\t" << "c_bulk\t" << "N\t" << "chi_s[kcal/mol]\t"
//			<< "chi_12[kcal/mol]\t" << endl;
//	cout << "A:\t" << A.fi_b << "\t" << A.N << "\t" << A.chi_s * R * temp
//			<< "\t\t" << A.chi_12 * R * temp << "\t" << endl;
//	cout << "B:\t" << B.fi_b << "\t" << B.N << "\t" << B.chi_s * R * temp
//			<< "\t\t" << B.chi_12 * R * temp << "\t" << endl;
//	cout << endl;
//
//	int fw = 8;
//	cout
//			<< "z   fi(A)   fi(B)   u(A)    u(B)    a(A)    a(B)    u_mix(A)u_mix(B)u_ads(A)u_ads(B)"
//			<< endl;
//	for (int z = 1; z < zmax; z++) {
//		cout << fixed << setprecision(3);
//		cout << setw(4) << left;
//		cout << z;
//		cout << setw(fw);
//		cout << A.fi_z[z];
//		cout << setw(fw);
//		cout << B.fi_z[z];
//		cout << setw(fw);
//		cout << A.u_total[z];
//		cout << setw(fw);
//		cout << B.u_total[z];
//		cout << setw(fw);
//		cout << A.alpha[z];
//		cout << setw(fw);
//		cout << B.alpha[z];
//		//cout<<A.u_adsorption[z]<<"\t";
//		cout << setw(fw);
//		cout << A.u_mixing[z];
//		cout << setw(fw);
//		cout << B.u_mixing[z];
//		cout << setw(fw);
//		cout << A.u_adsorption[z];
//		cout << setw(fw);
//		cout << B.u_adsorption[z] << endl;
//		//cout<<B.u_adsorption[z]<<"\t";
//
//	}
//
//	if (toFile == true) {
//		ofstream f;
//		f.open("data", ios::ate);
//		for (int z = 1; z < zmax; z++) {
//			f << fixed << setprecision(3);
//			f << setw(4) << left;
//			f << z;
//			f << setw(fw);
//			f << A.fi_z[z];
//			f << setw(fw);
//			f << B.fi_z[z];
//			f << setw(fw);
//			f << A.u_total[z];
//			f << setw(fw);
//			f << B.u_total[z];
//			f << setw(fw);
//			f << A.alpha[z];
//			f << setw(fw);
//			f << B.alpha[z];
//			//f<<A.u_adsorption[z]<<"\t";
//			f << setw(fw);
//			f << A.u_mixing[z];
//			f << setw(fw);
//			f << B.u_mixing[z];
//			f << setw(fw);
//			f << A.u_adsorption[z];
//			f << setw(fw);
//			f << B.u_adsorption[z] << endl;
//			//f<<B.u_adsorption[z]<<"\t";
//
//		}
//		f.close();
//	}
//
//}

void ScfCalc::printInfo3C(const Molecule& A, const Molecule& B,
		const Molecule& O, const double temp, bool verbose, bool toFile) {
	cout << fixed << setprecision(2);
	cout << endl << "Printing general info:" << endl;
	cout << "Species\t" << "c_bulk\t" << "N\t" << "chi_s[kcal/mol]\t"
			<< "chi_12[kcal/mol]" << "chi_O[kcal/mol]\t" << endl;
	cout << "A:\t" << A.fi_b << "\t" << A.N << "\t" << A.chi_s * R * temp
			<< "\t\t" << A.chi_12 * R * temp << "\t\t" << A.chi_O * R * temp
			<< "\t" << endl;
	cout << "B:\t" << B.fi_b << "\t" << B.N << "\t" << B.chi_s * R * temp
			<< "\t\t" << B.chi_12 * R * temp << "\t\t" << B.chi_O * R * temp
			<< "\t" << endl;
	cout << "Solvent\t" << "c_bulk\t" << "N\t" << "chi_s[kcal/mol]\t"
			<< "chi_12[kcal/mol]" << "chi_13[kcal/mol]\t" << endl;
	cout << "O:\t" << O.fi_b << "\t" << O.N << "\t" << O.chi_s * R * temp
			<< "\t\t" << O.chi_12 * R * temp << "\t\t" << O.chi_O * R * temp
			<< "\t" << endl;
	cout << endl;
	int fw = 8;
	if (verbose == true) {

		cout
				<< "z   fi(A)   fi(B)   fi(O)   u(A)    u(B)    u(O)    a(A)    a(B)    a(O)    u_mix(A)u_mix(B)u_mix(O)u_ads(A)u_ads(B)u_ads(O)"
				<< endl;
		for (int z = 1; z < zmax; z++) {
			cout << fixed << setprecision(3);
			cout << setw(4) << left;
			cout << z;
			cout << setw(fw);
			cout << A.fi_z[z];
			cout << setw(fw);
			cout << B.fi_z[z];
			cout << setw(fw);
			cout << O.fi_z[z];
			cout << setw(fw);
			cout << A.u_total[z];
			cout << setw(fw);
			cout << B.u_total[z];
			cout << setw(fw);
			cout << O.u_total[z];
			cout << setw(fw);
			cout << A.alpha[z];
			cout << setw(fw);
			cout << B.alpha[z];
			cout << setw(fw);
			cout << O.alpha[z];
			//cout<<A.u_adsorption[z]<<"\t";
			cout << setw(fw);
			cout << A.u_mixing[z];
			cout << setw(fw);
			cout << B.u_mixing[z];
			cout << setw(fw);
			cout << O.u_mixing[z];
			cout << setw(fw);
			cout << A.u_adsorption[z];
			cout << setw(fw);
			cout << B.u_adsorption[z];
			cout << setw(fw);
			cout << O.u_adsorption[z] << endl;
			//cout<<B.u_adsorption[z]<<"\t";

		}
	}
	if (toFile == true) {
		ofstream f;
		f.open("data", ios::ate);
		for (int z = 1; z < zmax; z++) {
			f << fixed << setprecision(3);
			f << setw(4) << left;
			f << z;
			f << setw(fw);
			f << A.fi_z[z];
			f << setw(fw);
			f << B.fi_z[z];
			f << setw(fw);
			f << O.fi_z[z];
			f << setw(fw);
			f << A.u_total[z];
			f << setw(fw);
			f << B.u_total[z];
			f << setw(fw);
			f << O.u_total[z];
			f << setw(fw);
			f << A.alpha[z];
			f << setw(fw);
			f << B.alpha[z];
			f << setw(fw);
			f << O.alpha[z];
			//f<<A.u_adsorption[z]<<"\t";
			f << setw(fw);
			f << A.u_mixing[z];
			f << setw(fw);
			f << B.u_mixing[z];
			f << setw(fw);
			f << O.u_mixing[z];
			f << setw(fw);
			f << A.u_adsorption[z];
			f << setw(fw);
			f << B.u_adsorption[z] << endl;
			f << setw(fw);
			f << O.u_adsorption[z] << endl;
			//f<<B.u_adsorption[z]<<"\t";
		}
		f.close();
	}

}

//void ScfCalc::printAdstats(const Molecule A, const Molecule B,
//		const double temp, bool toFile) {
//	cout << fixed << setprecision(1);
//	cout << endl << "Total amount A:\t" << setw(7) << right << A.total << " ("
//			<< (A.total / (A.total + B.total) * 100) << "%)" << endl;
//	cout << "Total amount B:\t" << setw(7) << right << B.total << " ("
//			<< (B.total / (A.total + B.total) * 100) << "%)" << endl;
//	cout << "Ads. amount A:\t" << setw(7) << right << A.adsorbed_amount << " ("
//			<< (A.adsorbed_amount / A.total * 100) << "%)" << endl;
//	cout << "Ads. amount B:\t" << setw(7) << right << B.adsorbed_amount << " ("
//			<< (B.adsorbed_amount / B.total * 100) << "%)" << endl;
//	cout << "Theta A:\t" << setw(4) << right << A.fi_z[1] / (A.fi_z[1]
//			+ B.fi_z[1]) * 100 << "%" << endl;
//	cout << "Theta B:\t" << setw(4) << right << B.fi_z[1] / (A.fi_z[1]
//			+ B.fi_z[1]) * 100 << "%" << endl << endl;
//	//	cout << fixed << setprecision(2);
//	//	cout <<"A fizfree[1]"<<A.fi_zfree[1];
//	//	cout <<" B fizfree[1]"<<B.fi_zfree[1];
//	//	cout <<" A fizads[1]"<<A.fi_zads[1];
//	//	cout <<" B fizads[1]"<<B.fi_zads[1];
//	//	cout <<" A fiz[1]"<<A.fi_z[1];
//	//	cout <<" B fiz[1]"<<B.fi_z[1];
//
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_ad", ios::app);
//		f << fixed << setprecision(4);
//		f << setw(6) << A.fi_b << setw(6) << setprecision(2)
//				<< A.adsorbed_amount << setw(6) << A.fi_z[1] / (A.fi_z[1]
//				+ B.fi_z[1]) << setw(8) << setprecision(4) << B.fi_b << setw(6)
//				<< setprecision(2) << B.adsorbed_amount << setw(6) << B.fi_z[1]
//				/ (A.fi_z[1] + B.fi_z[1]) << endl;
//		//f<<B.u_adsorption[z]<<"\t";
//
//
//		f.close();
//	}
//}

//void ScfCalc::printAdstats3C(const Molecule A, const Molecule B,
//		const Molecule O, const double temp, bool toFile) {
//	cout << fixed << setprecision(1);
//	cout << endl << "Total amount A:\t" << setw(7) << right << A.total << " ("
//			<< (A.total / (A.total + B.total + O.total) * 100) << "%)" << endl;
//	cout << "Total amount B:\t" << setw(7) << right << B.total << " ("
//			<< (B.total / (A.total + B.total + O.total) * 100) << "%)" << endl;
//	cout << "Total amount O:\t" << setw(7) << right << O.total << " ("
//			<< (O.total / (A.total + B.total + O.total) * 100) << "%)" << endl;
//	cout << "Ads. amount A:\t" << setw(7) << right << A.adsorbed_amount << " ("
//			<< (A.adsorbed_amount / A.total * 100) << "%)" << endl;
//	cout << "Ads. amount B:\t" << setw(7) << right << B.adsorbed_amount << " ("
//			<< (B.adsorbed_amount / B.total * 100) << "%)" << endl;
//	cout << "Ads. amount O:\t" << setw(7) << right << O.adsorbed_amount << " ("
//			<< (O.adsorbed_amount / O.total * 100) << "%)" << endl;
//	cout << "Theta A:\t" << setw(4) << right << A.fi_z[1] / (A.fi_z[1]
//			+ B.fi_z[1] + O.fi_z[1]) * 100 << "%" << endl;
//	cout << "Theta B:\t" << setw(4) << right << B.fi_z[1] / (A.fi_z[1]
//			+ B.fi_z[1] + O.fi_z[1]) * 100 << "%" << endl;
//	cout << "Theta O:\t" << setw(4) << right << O.fi_z[1] / (A.fi_z[1]
//			+ B.fi_z[1] + O.fi_z[1]) * 100 << "%" << endl << endl;
//	//	cout << fixed << setprecision(2);
//	//	cout <<"A fizfree[1]"<<A.fi_zfree[1];
//	//	cout <<" B fizfree[1]"<<B.fi_zfree[1];
//	//	cout <<" A fizads[1]"<<A.fi_zads[1];
//	//	cout <<" B fizads[1]"<<B.fi_zads[1];
//	//	cout <<" A fiz[1]"<<A.fi_z[1];
//	//	cout <<" B fiz[1]"<<B.fi_z[1];
//
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_ad", ios::app);
//		f << fixed << setprecision(4);
//		f << setw(6) << A.fi_b << setw(6) << setprecision(2)
//				<< A.adsorbed_amount << setw(6) << A.fi_z[1] / (A.fi_z[1]
//				+ B.fi_z[1]) << setw(8) << setprecision(4) << B.fi_b << setw(6)
//				<< setprecision(2) << B.adsorbed_amount << setw(6) << B.fi_z[1]
//				/ (A.fi_z[1] + B.fi_z[1]) << endl;
//		//f<<B.u_adsorption[z]<<"\t";
//
//
//		f.close();
//	}
//}

void ScfCalc::printAdstatsCopo(const Copolymer X, const Molecule O,
		const double temp, bool toFile) {
	cout << fixed << setprecision(2);
	cout << endl << "Total amount Copo:\t" << setw(7) << right << X.total
			<< " (" << (X.total / (X.total + O.total) * 100) << "%)" << endl;
	cout << "Excess amount Copo:\t" << setw(7) << right << (X.total
			- (zmax - 1) * X.fi_b) << endl;
	cout << "Total amount O:\t\t" << setw(7) << right << O.total << " ("
			<< (O.total / (X.total + O.total) * 100) << "%)" << endl << endl;
	cout << "Ads. amount Copo:\t" << setw(7) << right << X.adsorbed_amount
			<< " (" << (X.adsorbed_amount / X.total * 100) << "%)" << endl;
	cout << "Tail amount Copo:\t" << setw(7) << right << X.tail_amount << " ("
			<< (X.tail_amount / X.total * 100) << "%)" << endl;
	cout << "Loop amount Copo:\t" << setw(7) << right << X.loop_amount << " ("
			<< (X.loop_amount / X.total * 100) << "%)" << endl;
	cout << "Train amount Copo:\t" << setw(7) << right << X.fi_z[1] << " ("
			<< (X.fi_z[1] / X.total * 100) << "%)" << endl;
	cout << "Ads. amount O:\t\t" << setw(7) << right << O.adsorbed_amount
			<< " (" << (O.adsorbed_amount / O.total * 100) << "%)" << endl
			<< endl;

	cout << "Theta Copo:\t\t" << setw(4) << right << X.theta * 100 << "%"
			<< endl;
	cout << "Theta Segm. A:\t\t" << setw(4) << right << X.fi_A[1] / (X.fi_z[1]
			+ O.fi_z[1]) * 100 << "%" << endl;
	cout << "Theta Segm. B:\t\t" << setw(4) << right << X.fi_B[1] / (X.fi_z[1]
			+ O.fi_z[1]) * 100 << "%" << endl;
	cout << "Theta O:\t\t" << setw(4) << right << O.theta * 100 << "%" << endl
			<< endl;
	cout << "Volume/Surface ratio [V/S = zmax]: " << (zmax - 1) << endl;
	//	cout << "Fi_z[2]: " << X.fi_z[3] * 100 << " Fi_zfree[2]: " << X.fi_zfree[3]
	//			* 100 << " Fi_zads[1]: " << X.fi_zads[3] * 100 << " Fi_zads2[1]"
	//			<< X.fi_zads2[3] * 100 << endl;
	//	cout << fixed << setprecision(2);
	//	cout <<"A fizfree[1]"<<A.fi_zfree[1];
	//	cout <<" B fizfree[1]"<<B.fi_zfree[1];
	//	cout <<" A fizads[1]"<<A.fi_zads[1];
	//	cout <<" B fizads[1]"<<B.fi_zads[1];
	//	cout <<" A fiz[1]"<<A.fi_z[1];
	//	cout <<" B fiz[1]"<<B.fi_z[1];

	if (toFile == true) {
		ofstream f;
		f.open("data_ad", ios::app);
		//fi_b, ads_amount, theta, tail, loop, temp
		f << fixed << setprecision(4);
		f << setw(6) << X.fi_b;
		f << setw(6) << setprecision(2) << X.adsorbed_amount;
		//Theta
		f << setw(6) << X.theta;
		f << setw(8) << setprecision(4) << X.loop_amount;
		f << setw(8) << setprecision(4) << X.tail_amount;
		f << setw(8) << setprecision(2) << temp;
		f << " " << X.getSequence() << endl;
		//f<<B.u_adsorption[z]<<"\t";


		f.close();
	}
}

void ScfCalc::printAdstatsCopo_Ensemble(const std::vector<Copolymer> &lcoposet,
		const Molecule O, const double all, const double tamount,
		const double temp, bool toFile) {
	std::vector<Copolymer>::const_iterator it;
	cout << "\n###Adsorption Statistics###" << endl;
	cout << fixed << setprecision(2);
	for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
		cout << "\n>Copolymer " << it->ident;
		cout << endl << "Total amount:\t" << setw(7) << right << it->total
				<< " (" << (it->total / (all) * 100) << "%)" << endl;
		cout << "Excess amount:\t" << setw(7) << right << (it->total - (zmax
				- 1) * it->fi_b) << endl;
		cout << "Ads. amount:\t" << setw(7) << right << it->adsorbed_amount
				<< " (" << (it->adsorbed_amount / it->total * 100) << "%)"
				<< endl;
		cout << "Tail amount:\t" << setw(7) << right << it->tail_amount << " ("
				<< (it->tail_amount / it->total * 100) << "%)" << endl;
		cout << "Loop amount:\t" << setw(7) << right << it->loop_amount << " ("
				<< (it->loop_amount / it->total * 100) << "%)" << endl;
		cout << "Train amount:\t" << setw(7) << right << it->fi_z[1] << " ("
				<< (it->fi_z[1] / it->total * 100) << "%)" << endl;
		cout << "Surface cover.: " << setw(7) << right << it->theta * 100
				<< "%" << endl;
		cout << "Segm. A:\t" << setw(7) << right << it->fi_A[1] / (tamount)
				* 100 << "%" << endl;
		cout << "Segm. B:\t" << setw(7) << right << it->fi_B[1] / (tamount)
				* 100 << "%" << endl;
	}
	cout << "\n>Solvent(Homopolymer)" << endl;
	cout << "Total amount:\t" << setw(7) << right << O.total << " ("
			<< (O.total / (all) * 100) << "%)" << endl;
	cout << "Excess amount:\t" << setw(7) << right << (O.total - (zmax - 1)
			* O.fi_b) << endl;
	cout << "Ads. amount:\t" << setw(7) << right << O.adsorbed_amount << " ("
			<< (O.adsorbed_amount / O.total * 100) << "%)" << endl;
	cout << "Surface cover.: " << setw(7) << right << O.theta * 100 << "%"
			<< endl << endl;
	cout << "Volume/Surface ratio [V/S = zmax]: " << (zmax - 1) << endl;
	//	cout << "Fi_z[2]: " << X.fi_z[3] * 100 << " Fi_zfree[2]: " << X.fi_zfree[3]
	//			* 100 << " Fi_zads[1]: " << X.fi_zads[3] * 100 << " Fi_zads2[1]"
	//			<< X.fi_zads2[3] * 100 << endl;
	//	cout << fixed << setprecision(2);
	//	cout <<"A fizfree[1]"<<A.fi_zfree[1];
	//	cout <<" B fizfree[1]"<<B.fi_zfree[1];
	//	cout <<" A fizads[1]"<<A.fi_zads[1];
	//	cout <<" B fizads[1]"<<B.fi_zads[1];
	//	cout <<" A fiz[1]"<<A.fi_z[1];
	//	cout <<" B fiz[1]"<<B.fi_z[1];

	if (toFile == true) {
		ofstream f;
		f.open("data_ad", ios::app);
		for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
			// fi_b, adsorbed_amount, theta, loop, tail, sequence
			f << fixed << setprecision(4);
			f << setw(6) << it->fi_b;
			f << setw(6) << setprecision(2) << it->adsorbed_amount;
			f << setw(6) << it->theta;
			f << setw(8) << setprecision(4) << it->loop_amount;
			f << setw(8) << setprecision(4) << it->tail_amount;
			f << " " << it->getSequence() << " ";
		}
		f << setw(6) << O.fi_b;
		f << setw(6) << setprecision(2) << O.adsorbed_amount;
		//temp
		f << setw(8) << setprecision(2) << temp << endl;
		f.close();
	}
}

//void ScfCalc::printAdstatsCopo_Ensemble_old(const Copolymer X,
//		const Copolymer Y, const Molecule O, const double temp, bool toFile) {
//	cout << fixed << setprecision(2);
//	cout << endl << "Total amount Copo:\t" << setw(7) << right << X.total
//			<< " (" << (X.total / (X.total + Y.total + O.total) * 100) << "%)";
//	cout << setw(7) << right << Y.total << " (" << (Y.total / (X.total
//			+ Y.total + O.total) * 100) << "%)" << endl;
//	cout << "Excess amount Copo:\t" << setw(7) << right << (X.total
//			- (zmax - 1) * X.fi_b);
//	cout << setw(16) << right << (Y.total - (zmax - 1) * Y.fi_b) << endl;
//	cout << "Total amount O:\t\t" << setw(7) << right << O.total << " ("
//			<< (O.total / (X.total + O.total) * 100) << "%)" << endl;
//	cout << "Excess amount O:\t" << setw(7) << right << (O.total - (zmax - 1)
//			* O.fi_b) << endl << endl;
//	cout << "Ads. amount Copo:\t" << setw(7) << right << X.adsorbed_amount
//			<< " (" << (X.adsorbed_amount / X.total * 100) << "%)";
//	cout << setw(7) << right << Y.adsorbed_amount << " (" << (Y.adsorbed_amount
//			/ Y.total * 100) << "%)" << endl;
//	cout << "Tail amount Copo:\t" << setw(7) << right << X.tail_amount << " ("
//			<< (X.tail_amount / X.total * 100) << "%)";
//	cout << setw(7) << right << Y.tail_amount << " (" << (Y.tail_amount
//			/ Y.total * 100) << "%)" << endl;
//	cout << "Loop amount Copo:\t" << setw(7) << right << X.loop_amount << " ("
//			<< (X.loop_amount / X.total * 100) << "%)";
//	cout << setw(7) << right << Y.loop_amount << " (" << (Y.loop_amount
//			/ Y.total * 100) << "%)" << endl;
//	cout << "Train amount Copo:\t" << setw(7) << right << X.fi_z[1] << " ("
//			<< (X.fi_z[1] / X.total * 100) << "%)";
//	cout << setw(7) << right << Y.fi_z[1] << " ("
//			<< (Y.fi_z[1] / Y.total * 100) << "%)" << endl;
//	cout << "Ads. amount O:\t\t" << setw(7) << right << O.adsorbed_amount
//			<< " (" << (O.adsorbed_amount / O.total * 100) << "%)" << endl
//			<< endl;
//
//	cout << "Theta Copo:\t\t" << setw(7) << right << X.theta * 100 << "%"
//			<< setw(7) << right << Y.theta * 100 << "%" << endl;
//	cout << "Theta Segm. A:\t\t" << setw(7) << right << X.fi_A[1] / (X.fi_z[1]
//			+ Y.fi_z[1] + O.fi_z[1]) * 100 << "%";
//	cout << setw(7) << right << Y.fi_A[1] / (X.fi_z[1] + Y.fi_z[1] + O.fi_z[1])
//			* 100 << "%" << endl;
//	cout << "Theta Segm. B:\t\t" << setw(7) << right << X.fi_B[1] / (X.fi_z[1]
//			+ Y.fi_z[1] + O.fi_z[1]) * 100 << "%";
//	cout << setw(7) << right << Y.fi_B[1] / (X.fi_z[1] + Y.fi_z[1] + O.fi_z[1])
//			* 100 << "%" << endl;
//	cout << "Theta O:\t\t" << setw(7) << right << O.theta * 100 << "%" << endl
//			<< endl;
//	cout << "Volume/Surface ratio [V/S = zmax]: " << (zmax - 1) << endl;
//	//	cout << "Fi_z[2]: " << X.fi_z[3] * 100 << " Fi_zfree[2]: " << X.fi_zfree[3]
//	//			* 100 << " Fi_zads[1]: " << X.fi_zads[3] * 100 << " Fi_zads2[1]"
//	//			<< X.fi_zads2[3] * 100 << endl;
//	//	cout << fixed << setprecision(2);
//	//	cout <<"A fizfree[1]"<<A.fi_zfree[1];
//	//	cout <<" B fizfree[1]"<<B.fi_zfree[1];
//	//	cout <<" A fizads[1]"<<A.fi_zads[1];
//	//	cout <<" B fizads[1]"<<B.fi_zads[1];
//	//	cout <<" A fiz[1]"<<A.fi_z[1];
//	//	cout <<" B fiz[1]"<<B.fi_z[1];
//
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_ad", ios::app);
//		//fi_b, ads_amount, theta, tail, loop, temp
//		f << fixed << setprecision(4);
//		f << setw(6) << X.fi_b;
//		f << setw(6) << setprecision(2) << X.adsorbed_amount;
//		//Theta X
//		f << setw(6) << X.theta;
//		f << setw(8) << setprecision(4) << X.loop_amount;
//		f << setw(8) << setprecision(4) << X.tail_amount;
//		f << setw(8) << setprecision(2) << temp;
//		f << " " << X.getSequence();
//		f << setw(6) << Y.fi_b;
//		f << setw(6) << setprecision(2) << Y.adsorbed_amount;
//		//Theta Y
//		f << setw(6) << Y.theta;
//		f << setw(8) << setprecision(4) << Y.loop_amount;
//		f << setw(8) << setprecision(4) << Y.tail_amount;
//		f << setw(8) << setprecision(2) << temp;
//		f << " " << Y.getSequence() << endl;
//		//f<<B.u_adsorption[z]<<"\t";
//
//
//		f.close();
//	}
//}

//void ScfCalc::printfiz(const Molecule A, const Molecule B, bool toFile) {
//	int fw = 8;
//	for (int z = 1; z < zmax; z++) {
//		cout << fixed << setprecision(3);
//		cout << setw(4) << left;
//		cout << z;
//		cout << setw(fw);
//		cout << A.fi_z[z] + B.fi_z[z];
//		cout << setw(fw);
//		cout << A.fi_z[z];
//		cout << setw(fw);
//		cout << B.fi_z[z];
//		cout << setw(fw);
//		cout << A.fi_zfree[z];
//		cout << setw(fw);
//		cout << B.fi_zfree[z];
//		cout << setw(fw);
//		cout << A.fi_zads[z];
//		cout << setw(fw);
//		cout << B.fi_zads[z] << endl;
//		//cout<<B.u_adsorption[z]<<"\t";
//	}
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_fi", ios::ate);
//		for (int z = 1; z < zmax; z++) {
//			f << fixed << setprecision(3);
//			f << setw(4) << left;
//			f << z;
//			f << setw(fw);
//			f << A.fi_z[z] + B.fi_z[z];
//			f << setw(fw);
//			f << A.fi_z[z];
//			f << setw(fw);
//			f << B.fi_z[z];
//			f << setw(fw);
//			f << A.fi_zfree[z];
//			f << setw(fw);
//			f << B.fi_zfree[z];
//			f << setw(fw);
//			f << A.fi_zads[z];
//			f << setw(fw);
//			f << B.fi_zads[z] << endl;
//
//		}
//		f.close();
//	}
//
//}

//void ScfCalc::printfiz3C(const Molecule A, const Molecule B, const Molecule O,
//		bool toFile) {
//	int fw = 8;
//	for (int z = 1; z < zmax; z++) {
//		cout << fixed << setprecision(3);
//		cout << setw(4) << left;
//		cout << z;
//		cout << setw(fw);
//		cout << A.fi_z[z] + B.fi_z[z] + O.fi_z[z];
//		cout << setw(fw);
//		cout << A.fi_z[z];
//		cout << setw(fw);
//		cout << B.fi_z[z];
//		cout << setw(fw);
//		cout << O.fi_z[z];
//		cout << setw(fw);
//		cout << A.fi_zfree[z];
//		cout << setw(fw);
//		cout << B.fi_zfree[z];
//		cout << setw(fw);
//		cout << O.fi_zfree[z];
//		cout << setw(fw);
//		cout << A.fi_zads[z];
//		cout << setw(fw);
//		cout << B.fi_zads[z];
//		cout << setw(fw);
//		cout << O.fi_zads[z] << endl;
//		//cout<<B.u_adsorption[z]<<"\t";
//	}
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_fi", ios::ate);
//		for (int z = 1; z < zmax; z++) {
//			f << fixed << setprecision(3);
//			f << setw(4) << left;
//			f << z;
//			f << setw(fw);
//			f << A.fi_z[z] + B.fi_z[z] + O.fi_z[z];
//			f << setw(fw);
//			f << A.fi_z[z];
//			f << setw(fw);
//			f << B.fi_z[z];
//			f << setw(fw);
//			f << O.fi_z[z];
//			f << setw(fw);
//			f << A.fi_zfree[z];
//			f << setw(fw);
//			f << B.fi_zfree[z];
//			f << setw(fw);
//			f << O.fi_zfree[z];
//			f << setw(fw);
//			f << A.fi_zads[z];
//			f << setw(fw);
//			f << B.fi_zads[z];
//			f << setw(fw);
//			f << O.fi_zads[z] << endl;
//		}
//		f.close();
//	}
//
//}

//prints most important parameter fi againt z
void ScfCalc::printfizCopo(const Copolymer X, const Molecule O, bool toFile) {
	int fw = 8;
	for (int z = 1; z < zmax; z++) {
		cout << fixed << setprecision(3);
		cout << setw(fw) << left << (double) z * mlength;
		cout << setw(fw) << X.fi_z[z] + O.fi_z[z];
		cout << setw(fw) << X.fi_z[z];
		cout << setw(fw) << X.fi_A[z];
		cout << setw(fw) << X.fi_B[z];
		cout << setw(fw) << O.fi_z[z];
		cout << setw(fw) << X.fi_zfree[z];
		cout << setw(fw) << O.fi_zfree[z];
		cout << setw(fw) << X.fi_zads[z];
		cout << setw(fw) << O.fi_zads[z];
		cout << setw(fw) << X.fi_tail[z];
		cout << setw(fw) << X.fi_loop[z] << endl;
		//cout<<B.u_adsorption[z]<<"\t";
	}
	if (toFile == true) {
		ofstream f;
		f.open("data_fi", ios::ate);
		for (int z = 1; z < zmax; z++) {
			f << fixed << setprecision(3);
			f << setw(fw) << left << (double) z * mlength;
			f << setw(fw) << X.fi_z[z] + O.fi_z[z];
			f << setw(fw) << X.fi_z[z];
			f << setw(fw) << X.fi_A[z];
			f << setw(fw) << X.fi_B[z];
			f << setw(fw) << O.fi_z[z];
			f << setw(fw) << X.fi_zfree[z];
			f << setw(fw) << O.fi_zfree[z];
			f << setw(fw) << X.fi_zads[z];
			f << setw(fw) << O.fi_zads[z];
			f << setw(fw) << X.fi_tail[z];
			f << setw(fw) << X.fi_loop[z] << endl;
		}
		f.close();
	}

}

//prints most important parameter fi againt z
void ScfCalc::printfizCopo_Ensemble(const std::vector<Copolymer> &lcoposet,
		const Molecule O, bool toFile) {
	int fw = 8;
	double copo_total[zmax];
	std::vector<Copolymer>::const_iterator it;
	for (int z = 1; z < zmax; z++) {
		copo_total[z] = 0;
		//1
		cout << setw(fw) << left << ((double) z * mlength) - mlength;
		for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
			//(ident-1)*7+1;
			cout << fixed << setprecision(3);
			//cout << setw(fw) << X.fi_z[z] + Y.fi_z[z] + O.fi_z[z];
			//(ident-1)*7+2; etc
			cout << setw(fw) << it->fi_z[z];
			cout << setw(fw) << it->fi_A[z];
			cout << setw(fw) << it->fi_B[z];
			cout << setw(fw) << it->fi_zfree[z];
			cout << setw(fw) << it->fi_zads[z];
			cout << setw(fw) << it->fi_tail[z];
			cout << setw(fw) << it->fi_loop[z];
			copo_total[z] = copo_total[z] + it->fi_z[z];
		}
		cout << setw(fw) << O.fi_z[z];
		cout << setw(fw) << O.fi_zfree[z];
		cout << setw(fw) << O.fi_zads[z];
		cout << setw(fw) << copo_total[z] << endl;
	}

	if (toFile == true) {
		ofstream f;
		f.open("data_fi", ios::ate);
		for (int z = 1; z < zmax; z++) {
			f << fixed << setprecision(3);
			f << setw(fw) << left << ((double) z * mlength) - mlength;
			for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
				f << setw(fw) << it->fi_z[z];
				f << setw(fw) << it->fi_A[z];
				f << setw(fw) << it->fi_B[z];
				f << setw(fw) << it->fi_zfree[z];
				f << setw(fw) << it->fi_zads[z];
				f << setw(fw) << it->fi_tail[z];
				f << setw(fw) << it->fi_loop[z];
			}
			f << setw(fw) << O.fi_z[z];
			f << setw(fw) << O.fi_zfree[z];
			f << setw(fw) << O.fi_zads[z];
			f << setw(fw) << copo_total[z] << endl;
		}

		f.close();
	}

}

//prints most important parameter fi againt z
//use script pd_ens3d
void ScfCalc::printfizCopo_Ensemble3D(const std::vector<Copolymer> &lcoposet) {
	int fw = 8;
	std::vector<Copolymer>::const_iterator it;
	ofstream f;
	f.open("data_fi3D", ios::ate);
	f << fixed << setprecision(3);
	//dummy line
	for (int z = 1; z < zmax; z++) {
		f << setw(fw) << left << ((double) z * mlength) - mlength;
		f << setw(fw) << 0;
		//f << setw(fw) << it->fi_z[z]<<endl;
		f << setw(fw) << 0.0 << endl;
	}
	f << endl;
	for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
		for (int z = 1; z < zmax; z++) {
			f << setw(fw) << left << ((double) z * mlength) - mlength;
			//f << setw(fw) << it->N;
			//manual correction
			f << setw(fw) << (it->N/3);
			//f << setw(fw) << it->fi_z[z]<<endl;
			f << setw(fw) << it->fi_zads[z] << endl;
		}
		//sometimes useful for gnuplot
		f << endl;
	}

	f.close();

}
//prints most important parameter fi againt z
//void ScfCalc::printfizCopo_Ensemble_old(const Copolymer X, const Copolymer Y,
//		const Molecule O, bool toFile) {
//	int fw = 8;
//	for (int z = 1; z < zmax; z++) {
//		cout << fixed << setprecision(3);
//		cout << setw(fw) << left << (double) z * mlength;
//		cout << setw(fw) << X.fi_z[z] + Y.fi_z[z] + O.fi_z[z];
//		cout << setw(fw) << X.fi_z[z];
//		cout << setw(fw) << X.fi_A[z];
//		cout << setw(fw) << X.fi_B[z];
//		cout << setw(fw) << Y.fi_z[z];
//		cout << setw(fw) << Y.fi_A[z];
//		cout << setw(fw) << Y.fi_B[z];
//		cout << setw(fw) << O.fi_z[z];
//		cout << setw(fw) << X.fi_zfree[z];
//		cout << setw(fw) << Y.fi_zfree[z];
//		cout << setw(fw) << O.fi_zfree[z];
//		cout << setw(fw) << X.fi_zads[z];
//		cout << setw(fw) << Y.fi_zads[z];
//		cout << setw(fw) << O.fi_zads[z];
//		cout << setw(fw) << X.fi_tail[z];
//		cout << setw(fw) << Y.fi_tail[z];
//		cout << setw(fw) << X.fi_loop[z];
//		cout << setw(fw) << Y.fi_loop[z] << endl;
//		//cout<<B.u_adsorption[z]<<"\t";
//	}
//	if (toFile == true) {
//		ofstream f;
//		f.open("data_fi", ios::ate);
//		for (int z = 1; z < zmax; z++) {
//			f << fixed << setprecision(3);
//			f << setw(fw) << left << (double) z * mlength;
//			f << setw(fw) << X.fi_z[z] + Y.fi_z[z] + O.fi_z[z];
//			f << setw(fw) << X.fi_z[z];
//			f << setw(fw) << X.fi_A[z];
//			f << setw(fw) << X.fi_B[z];
//			//6
//			f << setw(fw) << Y.fi_z[z];
//			f << setw(fw) << Y.fi_A[z];
//			f << setw(fw) << Y.fi_B[z];
//			f << setw(fw) << O.fi_z[z];
//			//10
//			f << setw(fw) << X.fi_zfree[z];
//			f << setw(fw) << Y.fi_zfree[z];
//			f << setw(fw) << O.fi_zfree[z];
//			f << setw(fw) << X.fi_zads[z];
//			f << setw(fw) << Y.fi_zads[z];
//			//15
//			f << setw(fw) << O.fi_zads[z];
//			f << setw(fw) << X.fi_tail[z];
//			f << setw(fw) << Y.fi_tail[z];
//			f << setw(fw) << X.fi_loop[z];
//			f << setw(fw) << Y.fi_loop[z] << endl;
//		}
//		f.close();
//	}
//
//}

void ScfCalc::printfis(Copolymer X, bool toStdout, bool toFile) {
	//int fw = 8;
	if (toStdout == true) {
		for (int z = 1; z < zmax; z++) {
			for (int s = 1; s <= X.N; s++) {
				cout << fixed << setprecision(3);
				cout << setw(8) << left;
				cout << (double) z * mlength;
				cout << setw(4) << left;
				cout << s;
				cout << setw(4) << left;
				cout << X.getfi(z, s) << endl;
			}

			//cout<<B.u_adsorption[z]<<"\t";
		}
	}
	if (toFile == true) {
		ofstream f;
		f.open("data_fis", ios::ate);
		for (int z = 1; z < zmax; z++) {
			for (int s = 1; s <= X.N; s++) {
				f << fixed << setprecision(3);
				f << setw(8) << left;
				f << (double) z * mlength;
				f << setw(4) << left;
				f << s;
				f << setw(4) << left;
				f << X.getfi(z, s) << endl;
			}

		}
		f.close();

	}
}
//Saving final segment density of a calculation
void ScfCalc::savePotential(const Molecule A, const Molecule B,
		const Molecule O) {
	for (int z = 1; z < zmax; z++) {
		fi_zA[z] = A.fi_z[z];
		fi_zB[z] = B.fi_z[z];
		fi_zO[z] = O.fi_z[z];
		u_A[z] = A.u_total[z];
		u_B[z] = B.u_total[z];
		u_O[z] = O.u_total[z];
	}
	fi_exist = true;
}
//Reading priorily saved segment densities
void ScfCalc::setPotential(Copolymer& X, Molecule& A, Molecule &B, Molecule& O,
		double temp) {
	for (int z = 1; z < zmax; z++) {
		A.fi_z[z] = fi_zA[z];
		B.fi_z[z] = fi_zB[z];
		O.fi_z[z] = fi_zO[z];
		A.u_total[z] = u_A[z] / 2.0;
		B.u_total[z] = u_B[z] / 2.0;
		O.u_total[z] = u_O[z] / 2.0;
	}
	X.calcG(1, temp, A.u_total, B.u_total);
	X.calcG_N(1, temp, A.u_total, B.u_total);
	X.calcfi(temp, A.u_total, B.u_total);
	//Transfer fi_A from Copolymer to the monomers
	for (int z = 1; z < zmax; z++) {
		////		//Something is going wrong with G, could it be because of G_avg
		////		cout << "A.fi_z["<<z<<"] "<<A.fi_z[z];
		////		//with function!
		A.fi_z[z] = X.fi_A[z];
		B.fi_z[z] = X.fi_B[z];
		////		cout << " A.fi_z["<<z<<"] "<<A.fi_z[z]<<endl;
		////		//cout << " B.fi_z["<<z<<"] "<<(*B).fi_z[z]<<endl;
	}
	A.calcG(1, temp);
	A.calcG_N(1, temp);
	A.calcfi();
	B.calcG(1, temp);
	B.calcG_N(1, temp);
	B.calcfi();
	O.calcG(1, temp);
	O.calcG_N(1, temp);
	O.calcfi();

}
//creates new sequence
void ScfCalc::shuffleSequence(const string oldseq, string& newseq) {
	//making a temporary copy
	for (int s = 0; s < (int) oldseq.length(); s++) {
		newseq[s] = oldseq[s];
	}
	cout << "Shuffle sequence..." << endl;
	int newpos;
	for (int s = 0; s < (int) oldseq.length(); s++) {
		//Polseq[s] = lPolseq[s];
		//newseq[s]=oldseq[s];
		//newpos=(int) (myRan.doub()*oldseq.length());
		//newpos = (int) (gsl_ran_flat(r, 0, 1) * oldseq.length());
		double number = distribution(r);
		newpos = (int) (number * oldseq.length());
		//cout <<"newpos: "<<newpos <<" s: "<<s<<" Oldseq[s]:"<<oldseq[s]<<" Oldseq[newpos]:"<<oldseq[newpos]<<endl;
		//switching
		if (newseq[s] == 'A' && newseq[newpos] == 'B') {
			newseq[s] = 'B';
			newseq[newpos] = 'A';
			//cout <<"newpos: "<<newpos <<" s: "<<s<<" Newseq[s]:"<<newseq[s]<<" Newseq[newpos]:"<<newseq[newpos]<<endl;
		}
		if (newseq[s] == 'B' && newseq[newpos] == 'A') {
			newseq[s] = 'A';
			newseq[newpos] = 'B';
			//		cout <<"newpos: "<<newpos <<" s: "<<s<<" Newseq[s]:"<<newseq[s]<<" Newseq[newpos]:"<<newseq[newpos]<<endl;
		}
	}

}

//Shuffle whole blocks
void ScfCalc::shuffleBlocks(const string block1, int numb1,
		const string block2, int numb2, string& newseq) {
	int total_length, number_blocks, newpos;
	bool pos_exists = false;
	total_length = block1.length() * numb1 + block2.length() * numb2;
	number_blocks = numb1 + numb2;
	string tempseq[number_blocks];
	int random_order[number_blocks];
	cout << "Total length: " << total_length << " Number of Blocks: "
			<< number_blocks << endl;
	newseq.clear();

	//putting block into string array
	for (int i = 0; i < numb1; i++) {
		//cout << i<< endl;
		tempseq[i] = block1;
	}
	for (int i = numb1; i < number_blocks; i++) {
		tempseq[i] = block2;
		//cout << i<< endl;
	}

	//Determine random order
	for (int i = 0; i < number_blocks; i++) {
		//generates random integer in between 1 and number_blocks
		do {
			double number = distribution(r);
			newpos = (int) (number * number_blocks + 1);
			pos_exists = false;
			//cout << "newpos: " << newpos << " ";
			for (int j = 0; j <= i; j++) {
				if (newpos == random_order[j]) {
					//cout << " newpos already at position: " << j;
					pos_exists = true;
				}
			}
		} while (pos_exists == true);
		//newpos = (int) (gsl_ran_flat(r, 0, 1) * number_blocks + 1);
		//cout << " new newpos: " << newpos << endl;
		pos_exists = false;
		random_order[i] = newpos;
		//cout << " random_array:" << random_order[i] << endl;
		//creating random block sequence
		newseq.append(tempseq[newpos - 1]);
		pos_exists = false;
	}

	//	newseq.append(block1);
	//	newseq.append(block2);
	//	newseq.append(block1);
}

//Function prints timing
void ScfCalc::printTiming(timeval &start, timeval &end) {
	int difsec;
	int difusec;
	difsec = end.tv_sec - start.tv_sec;
	difusec = end.tv_usec - start.tv_usec;
	//cout << start.tv_sec << ':' << start.tv_usec << endl;
	//cout << end.tv_sec << ':' << end.tv_usec << endl;
	//cout <<"Simulation ended after: " << setprecision(6) << dif << " sec.";
	if (difusec < 0) {
		difsec--;
		difusec = 1000000 + difusec;
		//cout <<"modifying";
	}
	cout << "\nTIMING: ";
	cout << fixed << setprecision(2);
	cout << difsec << "." << difusec << " sec" << endl;
}
//reads distributions
void ScfCalc::readDist(std::vector<string>& seqset,
		std::vector<double>& concset) {
	ifstream myfile1;
	std::string temp, temp3;
	std::string line;
	char *temp2;
	double conc = 0.0;
	myfile1.open("distribution");
	cout << "Reading polymer distribution from file..." << endl;
	const boost::regex seq("([A|a|B|b]{1,})[[:blank:]]+([0-9]{1,}\\.?[0-9]*)");
	//const boost::regex short_seq("\\(([A|a|B|b]{1,})\\)([0-9]){1,}([[:blank:]]+[0-9]{1,}\\.?[0-9]*)");
	const boost::regex
			short_seq(
					"\\(([A|a|B|b]{1,})\\)([0-9]{1,})([[:blank:]]+[0-9]{1,}\\.?[0-9]*)");
	if (myfile1.is_open()) {
		while (!myfile1.eof()) {
			getline(myfile1, line);
			boost::smatch matches, matches_short;
			//substitue short notation

				if (boost::regex_search(line, matches_short, short_seq)) {
					cout << "Found short notation: " << matches_short[1].str()
							<< " " << matches_short[2].str() << " times"
							<< endl;
					int n = 0;
					temp = matches_short[2].str();
					temp2 = new char[temp.length() + 1];
					strcpy(temp2, temp.c_str());
					n = atoi(temp2);
					temp3.append("\\3");
					for (int i = 0; i < n; i++) {
						temp3.insert(0, "\\1");
					}
					//cout<<temp3;
					//line.clear();
					//line.append(temp3);
					//line.append(matches_short[3].str());
					const std::string long_notation(temp3);
					line = boost::regex_replace(line, short_seq, long_notation,
							boost::match_default | boost::format_sed);
					temp3.clear();
				}

			if (boost::regex_search(line, matches, seq)) {
				cout << "Sequence: " << matches[1].str();
				seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				conc = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << "\nvol.fraction: " << conc << endl;
			}

			//cout << line;
		}
		myfile1.close();
	} else {
		cout << "Unable to open file ";
	}
	myfile1.close();
	//exit(1);
}

//gets new concentration according to adsorption profile (negative), normalized!
void ScfCalc::getadsDist(const std::vector<Copolymer>& lcoposet,
		const std::vector<double>& concset, std::vector<double>& adset) {
	//Hallo
	int i = 0;
	double sum = 0.0;
	std::vector<Copolymer>::const_iterator it;
	//get total adsorbed amount
	for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
		sum = sum + it->adsorbed_amount;
	}
	cout << fixed << setprecision(3);
	for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
		cout << "P#" << it->ident << " initial_c " << concset[i]
				<< " adsorbed_c ";
		cout << it->adsorbed_amount / sum;
		adset.push_back(it->adsorbed_amount / sum);
		cout << "  delta " << (adset[i] - concset[i]) << endl;
		i++;
	}
	cout << "Total adsorbed amount: " << sum << endl;
}
//normalizes distribution
void ScfCalc::normalizeDist(std::vector<double>& n_set, double nconc) {
	//Info
	double sum = 0;
	//std::vector<double>::const_iterator it;
	for (int i = 0; i < (int) n_set.size(); i++) {
		sum = sum + n_set[i];
	}
	cout << "Total amount: " << sum << "\tNormalisation of distribution to "
			<< nconc << endl;
	;
	for (int i = 0; i < (int) n_set.size(); i++) {
		n_set[i] = n_set[i] / sum * nconc;
		//cout << n_set[i]<<endl;
	}
	//exit(1);
}
//Writes distribution to file
void ScfCalc::writeDist(const std::vector<Copolymer>& lcoposet,
		const std::vector<double>& concset) {
	ofstream f;
	f.open("dist_adsorbed", ios::ate);
	std::vector<Copolymer>::const_iterator it;
	int i = 0;
	for (it = lcoposet.begin(); it != lcoposet.end(); ++it) {
		// fi_b, adsorbed_amount, theta, loop, tail, sequence
		f << it->getSequence() << "\t";
		f << fixed << setprecision(4);
		f << setw(6) << concset[i] << endl;
		i++;
	}
	f.close();
}

void ScfCalc::readParams(double& chi_AB, double& chi_As, double& chi_Bs,
		double& chi_Asol, double& chi_Bsol, double& temperature, double& csolvent) {
	ifstream myfile1;
	std::string line, temp;
	char *temp2;
	double tmp = 0.0;
	myfile1.open("SETUP");
	cout << "Reading polymer SETUP file..." << endl;
	cout << fixed << setprecision(2);
	boost::regex reg_chi_AB("(chi_ab)[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex reg_chi_As(
			"(chi_as)[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex reg_chi_Bs(
			"(chi_bs)[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex reg_chi_Asol(
			"(chi_asol)[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex reg_chi_Bsol(
				"(chi_bsol)[[:blank:]]+(-?[0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex
			reg_temp("(temp|TEMP|Temp|T)[[:blank:]]+([0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	boost::regex reg_csolv(
			"(csolv)[[:blank:]]+([0-9]{0,}\\.?[0-9]*)", boost::regex::icase);
	if (myfile1.is_open()) {
		while (!myfile1.eof()) {
			getline(myfile1, line);
			boost::smatch matches;
			if (boost::regex_search(line, matches, reg_chi_AB)) {
				cout << "A-B interaction:  \t chi_AB=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << " [*RT]" << endl;
				chi_AB = tmp;
			}
			if (boost::regex_search(line, matches, reg_chi_As)) {
				cout << "A-surface interaction:\t chi_As=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << " [*RT]" << endl;
				chi_As = tmp;
			}
			if (boost::regex_search(line, matches, reg_chi_Bs)) {
				cout << "B-surface interaction:\t chi_Bs=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << " [*RT]" << endl;
				chi_Bs = tmp;
			}
			if (boost::regex_search(line, matches, reg_chi_Asol)) {
				cout << "A-solvent interaction:\t chi_Asol=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << " [*RT]" << endl;
				chi_Asol = tmp;
			}
			if (boost::regex_search(line, matches, reg_chi_Bsol)) {
							cout << "B-solvent interaction:\t chi_Bsol=";
							//seqset.push_back(matches[1].str());
							//converting to double
							temp = matches[2].str();
							temp2 = new char[temp.length() + 1];
							strcpy(temp2, temp.c_str());
							tmp = atof(temp2);
							//cout <<"1: "<< matches[0].str();
							//concset.push_back(atof(temp2));
							//cout <<" 3: "<< matches[2].str();
							cout << tmp << " [*RT]" << endl;
							chi_Bsol = tmp;
						}
			if (boost::regex_search(line, matches, reg_csolv)) {
				cout << "Concentration solvent:\t csolv=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << endl;
				csolvent = tmp;
			}
			if (boost::regex_search(line, matches, reg_temp)) {
				cout << "Temperature=";
				//seqset.push_back(matches[1].str());
				//converting to double
				temp = matches[2].str();
				temp2 = new char[temp.length() + 1];
				strcpy(temp2, temp.c_str());
				tmp = atof(temp2);
				//cout <<"1: "<< matches[0].str();
				//concset.push_back(atof(temp2));
				//cout <<" 3: "<< matches[2].str();
				cout << tmp << " K" << endl;
				temperature = tmp;
			}

		}
		myfile1.close();
	} else {
		cout << "Unable to open file ";
	}
	cout << endl;
	//exit(1);
}

//Destructor
ScfCalc::~ScfCalc() {
	//gsl_rng_free(r);

	// TODO Auto-generated destructor stub
}

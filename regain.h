/*
 * =====================================================================================
 *
 *       Filename:  regain.h
 *
 *    Description:  Regression GAIN calculation.  Uses linear or logistic regression in 
 *    				calculating pairwise interaction between SNPs.
 *
 *        Created:  02/02/2012
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu
 *
 * =====================================================================================
 */

#ifndef __REGAIN_H__
#define __REGAIN_H__

#include "plink/zfstream.h"

using namespace std;

// type for storing p-value and matrix position (row,col) of
// reGAIN interaction terms
typedef pair< double, pair<int, int> > mat_el;

class Regain {
	public:
		Regain(bool compr, double sifthr, bool integrative, bool fdrpr = false);
		~Regain();
		virtual void run();
		void mainEffect(int e1, bool numeric);
		void addCovariates(Model &m);
		void interactionEffect(int e1, bool numeric1, int e2, bool numeric2);
		virtual void writeRegain(bool fdrprune = false);
		virtual void writePvals();
		virtual void fdrPrune(double fdr);
		void writeRcomm(double T, double fdr);
		static bool mecomp (const mat_el &l, const mat_el &r);

	private:
		// integrative regain mode
		bool intregain;
		// num attributes (SNPs + numeric for integrative, SNPs for normal regain)
		int numattr;
		// use zlib compression?
		bool compressed;
		// apply FDR pruning matrix?
		bool fdrprune;
		// SIF interaction threshold
		double sif_thresh;
		// Output matrix files
		ZOutput REGAIN_MATRIX;
		ZOutput REGAIN_PMATRIX;
		// in memory arrays
		double** gainMatrix;
		double** gainPMatrix;
		// addition output files
		ofstream MEBETAS;
		ofstream BETAS;
		ofstream SIF;
		// collection of all interaction terms as mat_el types
		vector<mat_el> gainPint;
};
#endif

/*
 * =====================================================================================
 *
 *       Filename:  regain.cpp
 *
 *    Description:  Regression GAIN calculation
 *
 *        Created:  06/20/2011
 *
 *        Author:  Nick Davis, nick-davis@utulsa.edu
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "plink/plink.h"
#include "plink/options.h"
#include "plink/model.h"
#include "plink/logistic.h"
#include "plink/linear.h"
#include "plink/stats.h"
#include "plink/helper.h"

#include "regain.h"

// Plink object
extern Plink* PP;

// constructor
Regain::Regain(bool compr, double sifthr, bool integrative, bool fdrpr) {
	// set class vars to passed args
	compressed = compr;
	sif_thresh = sifthr;
	intregain = integrative;
	fdrprune = fdrpr;

	// set integrative/normal regain vars
	// additional ext for integrative
	string ext = intregain ? ".int" : "";
	// header in betas files
	string hdr = intregain ? "attr" : "SNP";
	// total number of attributes
	numattr = intregain ? PP->nl_all + par::nlist_number : PP->nl_all;

	// initialize matrices and open output files
	string beta_f = par::output_file_name + ext + ".betas";
	string mebeta_f = par::output_file_name + ext + ".mebetas";
	BETAS.open(beta_f.c_str(), ios::out);
	MEBETAS.open(mebeta_f.c_str(), ios::out);
	cout << "Writing epistasis pairwise beta values to [ " << beta_f << " ]" << endl;
	cout << "Writing epistasis main effect beta values to [ " << mebeta_f << " ]" << endl;
	BETAS.precision(6);
	MEBETAS.precision(6);
	// print header
	BETAS << hdr << "1\t" << hdr << "2\tB_0\tB_1\tB_1 P-VAL\tB_2\tB_2 P-VAL";
	if(par::covar_file) {
		for (int i = 0; i < par::clist_number; i++) {
				BETAS << "\t" << PP->clistname[i] << "\t" << PP->clistname[i] << " P-VAL";
		}
	}
	BETAS << "\tB_3\tB_3 P-VAL" << endl;

	MEBETAS << hdr << "\tB_0\tB_1\tB_1 P-VAL";
	if(par::covar_file) {
		for (int i = 0; i < par::clist_number; i++) {
				MEBETAS << "\t" << PP->clistname[i] << "\t" << PP->clistname[i] << " P-VAL";
		}
	}
	MEBETAS << endl;

	string sif_f = par::output_file_name + ext + ".sif";
	SIF.open(sif_f.c_str(), ios::out);
	cout << "Writing SIF file to [ " << sif_f << " ]" << endl;
	SIF.precision(6);
	
	gainMatrix = new double*[numattr];
	gainPMatrix = new double*[numattr];
	// allocate reGAIN matrix
	for(int i=0; i < numattr; ++i) {
		gainMatrix[i] = new double[numattr];
	}
	// allocate reGAIN p-value matrix
	for(int i=0; i < numattr; ++i) {
		gainPMatrix[i] = new double[numattr];
	}
}

// construct regression models for each pair of SNPs
void Regain::run() {
  // Count how many items in the SET1
  int e1, e2;
#ifdef _OPENMP
  // OpenMP parallelization of this outer loop
  int numThreads = omp_get_num_threads();
  int numProcs = omp_get_num_procs();
  cout << "\t\t" << numThreads << " OpenMP threads available" << endl;
  cout << "\t\t" << numProcs << " OpenMP processors available" << endl;
  // omp_set_num_threads(numProcs);
  // printLOG("OpenMP number of threads set to    " + dbl2str((double) numThreads) + "\n");
  // printLOG("OpenMP number of processors set to " + dbl2str((double) numProcs) + "\n");
  // omp_set_nested(1);
  #pragma omp parallel for schedule(dynamic, 1) private(e1, e2)
#endif
  for(e1 = 0; e1 < numattr; e1++) {
//	cout << "Peforming tests of epistasis: group "
//			<< ++epcc << " of " << epc << "        \r";
//	cout.flush();

      for(e2 = 0; e2 < numattr; e2++) {
        // We've already performed this test, since the matrix is symmetric
        if(e1 > e2) continue;

		// main effect of SNP 1 - for I_2 diagonal of the GAIN matrix
        if(e1 == e2)
#ifdef _OPENMP
#pragma omp critical
#endif
			mainEffect(e1, e1 >= PP->nl_all);

		else interactionEffect(e1, e1 >= PP->nl_all, e2, e2 >= PP->nl_all);

      }

   } // Next pair of SNPs
}

// corresponds to I_2 in GAIN
void Regain::mainEffect(int e1, bool numeric){
	Model *lm_main_effect;

	if(par::bt) {
	  LogisticModel * m = new LogisticModel(PP);
	  lm_main_effect = m;
	} else {
	  LinearModel * m = new LinearModel(PP);
	  lm_main_effect = m;
	}

	// Set missing data
	lm_main_effect->setMissing();

	string label = numeric ? "NUM" : "ADD";

	// Main effect of SNP/numeric attribute
	if (numeric) lm_main_effect->addNumeric(e1 - PP->nl_all);
	else lm_main_effect->addAdditiveSNP(e1);
	lm_main_effect->label.push_back(label);

	// handle covariates for reGAIN
	if(par::covar_file) {
		addCovariates(*lm_main_effect);
	}

	// Build design matrix
	lm_main_effect->buildDesignMatrix();

	// Fit linear model
	lm_main_effect->fitLM();

	// Did model fit okay?
	lm_main_effect->validParameters();

	// Obtain estimates and statistic
	int tp = 1; // always use first coefficient after intercept as main effect term
	lm_main_effect->testParameter = tp; // single variable main effect
	vector_t b_main_effect = lm_main_effect->getCoefs();
	vector_t b_p_values = lm_main_effect->getPVals();

	// set main effect (diagonal) beta coefficient and corresponding
	// p-value
	gainMatrix[e1][e1] = b_main_effect[tp];
	gainPMatrix[e1][e1] = b_p_values[tp - 1]; // p-values don't include intercept term

	// update main effect betas file
	if (numeric) MEBETAS << PP->nlistname[e1 - PP->nl_all];
	else MEBETAS << PP->locus[e1]->name;
	for(unsigned int i=0; i < b_main_effect.size(); ++i) {
		if (i == 0) // B0 coefficient doesn't have pval
			MEBETAS << "\t" << b_main_effect[i];
		// adjust pvals index since there's no B0 pval
		else
			MEBETAS << "\t" << b_main_effect[i] << "\t" << b_p_values[i - 1];

	}
	MEBETAS << endl;
	
	delete lm_main_effect;
}

// handle multiple covariates
void Regain::addCovariates(Model &m){
	for ( int i = 0; i < par::clist_number; i++ ) {
		// add covariate to the model
		m.addCovariate(i);
		m.label.push_back(PP->clistname[i]);
	}
}

// corresponds to I_3 in GAIN
// make full symmetric matrix
void Regain::interactionEffect(int e1, bool numeric1, int e2, bool numeric2) {
	///////////////////////////////////////////////
	// Logistic or linear regression epistasis test
	Model * lm;

	if(par::bt) {
	LogisticModel * m = new LogisticModel(PP);
	lm = m;
	} else {
	LinearModel * m = new LinearModel(PP);
	lm = m;
	}

	// Set missing data
	lm->setMissing();

	string label1 = numeric1 ? "NUM1" : "ADD1";
	string label2 = numeric2 ? "NUM2" : "ADD2";

	// Main effect of SNP/numeric attribute 1
	if (numeric1) lm->addNumeric(e1 - PP->nl_all);
	else lm->addAdditiveSNP(e1);
	lm->label.push_back(label1);

	// Main effect of SNP/numeric attribute 2
	if (numeric2) lm->addNumeric(e2 - PP->nl_all);
	else lm->addAdditiveSNP(e2);
	lm->label.push_back(label2);

	// handle covariates for reGAIN
	if(par::covar_file) {
		addCovariates(*lm);
	}

	// Epistasis
	lm->addInteraction(1, 2);
	lm->label.push_back("EPI");

	// Build design matrix
	lm->buildDesignMatrix();

	// Prune out any remaining missing individuals
	// No longer needed
	//		   lm->pruneY();

	// Fit linear model
	lm->fitLM();

	// Did model fit okay?
	lm->validParameters();

	vector_t b;
	// interaction
	int tp = 3;

	// add # covars to test param to get interaction param
	if (par::covar_file){
		tp += par::clist_number;
	}
	lm->testParameter = tp; // interaction
#ifdef _OPENMP
#pragma omp critical
{
#endif
	b = lm->getCoefs();

	//printLOG("The model p-value (x^2) = " + dbl2str(pvalue) + "\n");
	
		vector_t beta_p = lm->getPVals();
	gainMatrix[e1][e2] = b[b.size() - 1];
	gainMatrix[e2][e1] = b[b.size() - 1];
	gainPMatrix[e1][e2] = beta_p[beta_p.size() - 1];
	gainPMatrix[e2][e1] = beta_p[beta_p.size() - 1];

	// store p-value along with (e1, e2) location of
	// item.  This is used later for FDR pruning
	if(fdrprune){
		pair<int, int> p = make_pair(e1,e2);
		mat_el Pint = make_pair(beta_p[beta_p.size() - 1], p);
		gainPint.push_back(Pint);
	}

	// update BETAS and SIF files

	// numeric attributes or SNPs
	if (numeric1) BETAS << PP->nlistname[e1 - PP->nl_all] << "\t";
	else BETAS << PP->locus[e1]->name << "\t";
	if (numeric2) BETAS << PP->nlistname[e2 - PP->nl_all];
	else BETAS << PP->locus[e2]->name;

	for(unsigned int i=0; i < b.size(); ++i) {
		if (i == 0) // B0 coefficient doesn't have pval
			BETAS << "\t" << b[i];
		// adjust pvals index since there's no B0 pval
		else
			BETAS << "\t" << b[i] << "\t" << beta_p[i - 1];

	}
	BETAS << endl;

	// add to SIF if interaction >= threshold
	if (abs(b[b.size() - 1]) >= sif_thresh) {
		if (numeric1)
			SIF << PP->nlistname[e1 - PP->nl_all] << "\t" << abs(b[b.size() - 1])  << "\t";
		else SIF << PP->locus[e1]->name << "\t" << abs(b[b.size() - 1])  << "\t";
		if (numeric2)
			SIF << PP->nlistname[e2 - PP->nl_all] << endl;
		else SIF << PP->locus[e2]->name << endl;
	}
#ifdef _OPENMP
	}
#endif
	// Clean up
	delete lm;
}

// write contents of reGAIN matrix to file.  If fdr is true, the filename is .pruned.regain, and
// the log output reflects this.
void Regain::writeRegain(bool fdrprune){
	// write the reGAIN matrix to file named <dataset>.regain
	string regain_matrix_f = par::output_file_name;
	// additional ext for integrative
	string ext = intregain ? ".int" : "";
	string tail = compressed ? ".gz" : "";
	if (fdrprune) {
		// close stream from previous write
		REGAIN_MATRIX.close();
//		REGAIN_MATRIX.clear();
		regain_matrix_f += ext + ".pruned.regain" + tail;
		cout << "Writing FDR-pruned epistasis REGAIN matrix [ " << regain_matrix_f << " ]" << endl;
	}
	else {
		regain_matrix_f += ext + ".regain" + tail;
		cout << "Writing epistasis REGAIN matrix [ " << regain_matrix_f << " ]" << endl;
	}
	REGAIN_MATRIX.open(regain_matrix_f.c_str(), compressed);
	// write SNP column names
	for(int cn=0; cn < PP->nl_all; ++cn) {
		if(cn) {
			REGAIN_MATRIX << "\t" << PP->locus[cn]->name;
		}
		else {
			REGAIN_MATRIX << PP->locus[cn]->name;
		}
	}
	// write numeric attribute column names
	for(int cn=0; cn < par::nlist_number; ++cn)
			REGAIN_MATRIX << "\t" << PP->nlistname[cn];
	REGAIN_MATRIX << "\n";
	// write matrix entries
	for(int i=0; i < numattr; ++i) {
		for(int j=i; j < numattr; ++j) {
			// use absolute value (magnitude) of betas 
			gainMatrix[i][j] = abs(gainMatrix[i][j]);	
			if (j == i) {// fill in symmetric entries, replacing j < i with tabs
				string tabs = "";
				for (int k = 0; k < j; k++)
					tabs += "\t";
				REGAIN_MATRIX << tabs << dbl2str_fixed(gainMatrix[i][j], 6);
			}
			else {
				REGAIN_MATRIX << "\t" << dbl2str_fixed(gainMatrix[i][j], 6);
			}
		}
		REGAIN_MATRIX << "\n";
	}
}

// write contents of reGAIN p-values matrix to file.
void Regain::writePvals(){
	// write the beta p-values to file named <dataset>.pvals.regain
	string REGAIN_PMATRIX_f = par::output_file_name;
	// additional ext for integrative
	string ext = intregain ? ".int" : "";
	string tail = compressed ? ".gz" : "";
	REGAIN_PMATRIX_f += ext + ".pvals.regain" + tail;
	REGAIN_PMATRIX.open(REGAIN_PMATRIX_f.c_str(), compressed);
	cout << "Writing epistasis REGAIN p-value matrix [ " << REGAIN_PMATRIX_f << " ]" << endl;
	// write SNP column names
	for(int cn=0; cn < PP->nl_all; ++cn) {
		if(cn) {
			REGAIN_PMATRIX << "\t" << PP->locus[cn]->name;
		}
		else {
			REGAIN_PMATRIX << PP->locus[cn]->name;
		}
	}
	// write numeric attribute column names
	for(int cn=0; cn < par::nlist_number; ++cn)
			REGAIN_PMATRIX << "\t" << PP->nlistname[cn];
	REGAIN_PMATRIX << "\n";
	// write matrix entries
	for(int i=0; i < numattr; ++i) {
		for(int j=i; j < numattr; ++j) {
			if (j == i) {// fill in symmetric entries, replacing j < i with tabs
				string tabs = "";
				for (int k = 0; k < j; k++)
					tabs += "\t";
				REGAIN_PMATRIX << tabs << dbl2str_fixed(gainPMatrix[i][j], 6);
			}
			else {
				REGAIN_PMATRIX << "\t" << dbl2str_fixed(gainPMatrix[i][j], 6);
			}
		}
		REGAIN_PMATRIX << "\n";
	}
}

// Benjamini Hochberg FDR pruning - removes interaction 
// terms from reGAIN matrix based on BH FDR threshold
// code based on method described in All of Statistics p. 167
void Regain::fdrPrune(double fdr){
	cout << "Calculating Benjamini Hochberg FDR for pruning" << endl;
	int m = gainPint.size();
	// sort gain interaction mal_el type by p-value, maintaining
	// gainPMatrix location (row, col) with sorted values
	sort(gainPint.begin(), gainPint.end(), Regain::mecomp);
	
	// use rough FDR (RFDR) to estimate alpha based on input FDR
	double alpha =  2 * m * fdr / (m + 1);
	int R = -1;
	// BH method
	for(int i = 0; i < m; i++) {
		double l = (i + 1) * alpha / (double) m;
		// test whether current p-value < current l
		if (gainPint[i].first < l){
			R = i;
		}
		else break;
	}

	// BH threshold condition not met with any p-values, so exit
	if (R == -1){
		cout << "No p-value meets BH threshold criteria, so nothing pruned" << endl;
		return;
	}

	// BH rejection threshold
	double T = gainPint[R].first;
	cout << "BH rejection threshold: T = " + dbl2str(T) << ", R = " + int2str(R) + "" << endl;
	cout << "Pruning reGAIN interaction terms with p-values <= T (" + dbl2str(T) + ")" << endl;

	// now prune (set to 0.0) all values at or below R index
	for (int i = 0; i <= R; i++) {
		pair<int,int> p = gainPint[i].second;
		// symmetric matrix, so set [e1][e2] and [e2][e1]
		gainMatrix[p.first][p.second] = 0.0;
		gainMatrix[p.second][p.first] = 0.0;
	}
	cout << "Pruned " + int2str(R + 1) + " values from reGAIN interaction terms" << endl;
	// use threshold to write R commands to generate FDR plot 
	writeRcomm(T, fdr);
}

// writing R commands to plot FDR graph
void Regain::writeRcomm(double T, double fdr){
	ofstream RCOMM;
	RCOMM.precision(6);
	string fdr_r_file = par::output_file_name + ".R";
	string betas_file = par::output_file_name + ".betas";
	cout << "Writing R commands to generate FDR plot [" << fdr_r_file << "]" << endl;

	RCOMM.open(fdr_r_file.c_str(), ios::out);
	RCOMM << "fdrvars <- read.delim(\"" << betas_file << "\")" << endl;
	RCOMM << "betas <- fdrvars$B_3" << endl;
	RCOMM << "pvals <- fdrvars$B_3.P.VAL" << endl;
	RCOMM << "T <- " << T << endl;
	RCOMM << "fdr <- " << fdr << endl;
	RCOMM << "plot(betas, -log10(pvals), type=\"n\")" << endl;
	RCOMM << "abline(h=-log10(T), col=\"green4\" lwd=3)" << endl;
	RCOMM << "accept <- which(-log10(pvals) > -log10(T))" << endl;
	RCOMM << "reject <- which(-log10(pvals) <= -log10(T))" << endl;
	RCOMM << "prnidx <- fdr * length(betas[accept])" << endl;
	RCOMM << "srtaccbetas <- sort(betas[accept])" << endl;
	RCOMM << "prnval <- srtaccbetas[prnidx]" << endl;
	RCOMM << "if(prnidx%%1!=0){" << endl;
	RCOMM << "prnval <- (srtaccbetas[floor(prnidx)] + srtaccbetas[ceiling(prnidx)]) / 2" << endl;
	RCOMM << "}" << endl;
	RCOMM << "prunex <- which(betas <= prnval)" << endl;
	RCOMM << "pruney <- which(-log10(pvals) > -log10(T))" << endl;
	RCOMM << "prune <- intersect(prunex, pruney)" << endl;
	RCOMM << "points(betas[accept], -log10(pvals[accept]), bg=\"green4\", pch=21)" << endl;
	RCOMM << "points(betas[reject], -log10(pvals[reject]), bg=\"blue\", pch=21)" << endl;
	RCOMM << "points(betas[prune], -log10(pvals[prune]), bg=\"red\", pch=21)" << endl;
	RCOMM << "abline(v=prnval), col=\"red\", lwd=3)" << endl;
	RCOMM << "title(\"Scatter plot of -log10 transformed p-values vs. regression betas\")" << endl;
	RCOMM << "legend(\"topleft\", inset=.05, title=\"Type\", c(\"Accepted\", \"Rejected\", \"Pruned\"), pch=c(21,21, 21), pt.bg=c(\"green4\", \"blue\", \"red\"))" << endl;
	RCOMM.close();
}

// comparison fnc for mat_el types
bool Regain::mecomp ( const mat_el &l, const mat_el &r) {
        return l.first < r.first;
}

// free memory and close file streams associated with reGAIN
Regain::~Regain(){
	// close BETAS and SIF ofstreams
	BETAS.close();
	SIF.close();

	// reGAIN matrix
	REGAIN_MATRIX.close();  // free gain matrix memory

	// free gain matrix memory
	for(int i=0; i < numattr; ++i) {
		delete [] gainMatrix[i];
	}
	delete [] gainMatrix;

	// reGAIN p-value matrix
	REGAIN_PMATRIX.close();  // free gain matrix memory

	// free gain matrix memory
	for(int i=0; i < numattr; ++i) {
		delete [] gainPMatrix[i];
	}
	delete [] gainPMatrix;
}

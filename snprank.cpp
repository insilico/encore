/*
 * =====================================================================================
 *
 *       Filename:  snprank.cpp
 *
 *    Description:  SNPrank - single nucleotide polymorphism (SNP) ranking algorithm
 *
 *					Uses a GAIN file, together with a damping factor gamma, 
 *					(default is .85), to compute SNPrank scores.
 *					Prints a series of rows containing the SNP name, SNPrank score,
 *					information gain, sorted in descending order by SNPrank.
 *
 *        Created:  08/11/2011
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu 
 *
 * =====================================================================================
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "plink/helper.h"

#include "snprank.h"

using namespace boost;
SNPrank::SNPrank() {}

SNPrank::SNPrank(string file) {
	readFile(file);
}

SNPrank::~SNPrank() {}

void SNPrank::readFile(string filename) {
	ZInput gainf(filename, compressed(filename));
	string line;
	string delimiter = "";
	// read first line (header)
	line = gainf.readLine();
	// set delimiter based on discovery of , or whitespace
	if (line.find(',') != string::npos){
		delimiter = ",";
	}

	else {
		delimiter = "\t ";
	}

	// tokenize header
	vector<string> tokHeader;
	split(tokHeader, line, is_any_of(delimiter));
	// create header from first line
	setHeader(tokHeader);	
	
	// initialize dimensions of data matrix
	size_t dim = tokHeader.size();
	data.set_size(dim, dim);

	// read numeric data
	int row = 0;

	for (row = 0; row < dim; row++) {
		line = gainf.readLine();
		vector<string> tokLine;
		//trim(line);
		split(tokLine, line, is_any_of(delimiter));
		for (size_t col = 0; col < dim; col++) {
			trim(tokLine[col]);
			tokLine[col] = ((tokLine[col] == "") ? "0" : tokLine[col]);
			data (row, col) = lexical_cast<double>(tokLine[col]);
		}
	}
}

void SNPrank::snprank(vector<string> names, mat G, double gamma, string outFile) {
	// ensure G is symmetric, reflect upper triangular to lower triangular
	mat symmG = symmatu(G);
	// G diagonal is information gain
	mat Gdiag = symmG.diag();

	double Gtrace = trace(symmG);

	//vector of column sums of G
	rowvec colsum = sum(symmG);

	// find the indices of non-zero elements
	uvec colsum_nzidx = find(colsum);
	
	// sparse matrix where the nonzero indices are filled with 1/colsum
	mat D = zeros<mat>(symmG.n_rows, symmG.n_cols);
	mat fillD = ones(colsum_nzidx.n_elem, 1) / colsum.elem(colsum_nzidx);

	for (int i = 0; i < fillD.size(); i++){
		D(colsum_nzidx[i], colsum_nzidx[i]) = fillD[i];
	}		


	// non-zero elements of colsum/d_j have (1 - gamma) in the numerator of 
	// the second term (Eq. 5 from SNPrank paper)
	mat T_nz = ones(1, symmG.n_cols);
	T_nz.elem(colsum_nzidx) = (1 - gamma) * ones<vec>(colsum_nzidx.n_elem);

	// Compute T, Markov chain transition matrix
	mat T = (symmG * D * gamma) + (Gdiag * T_nz) / Gtrace;

	// initial arbitrary vector 
	vec r(symmG.n_rows);
	r.fill(1.0 / symmG.n_rows);
	
	double threshold = 1.0E-4;
	double lambda = 0.0;
	bool converged = false;
	vec r_old = r;

	// if the absolute value of the difference between old and current r 
	// vector is < threshold, we have converged
	while(!converged) {
		r_old = r;
		r = T * r;
		
		// sum of r elements
		lambda = sum(r);
		
		// normalize eigenvector r so sum(r) == 1
		r = r / lambda;
		
		// check convergence, ensure all elements of r - r_old < threshold
		if (min((uvec)(abs(r - r_old) < threshold)) == 1){
			converged = true;
		}
		
		else converged = false;
	}

	// r indices sorted in descending order 
	uvec r_indices = sort_index(r, 1);

	// output r (SNPrank rankings) to file, truncating to 6 decimal places
	ofstream snpout(outFile.c_str());	
	snpout << "SNP\tSNPrank\tIG" << endl;
	int index = 0;
	
	for (int i = 0; i < r.n_elem; i++) {
		index = r_indices[i];
		snpout << names[index] << "\t" << fixed << scientific << r[index] << "\t" << symmG(index, index)
			<< endl;
	}

	snpout.close();
}

vector<string> SNPrank::getHeader() {
	return header;
}

void SNPrank::setHeader(vector<string> hdr) {
	header = hdr;
}

mat SNPrank::getData() {
	return data;
}

void SNPrank::setData(mat dat) {
	data = dat;
}

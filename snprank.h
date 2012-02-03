/*
 * =====================================================================================
 *
 *       Filename:  snprank.h
 *
 *    Description:  SNPrank - single nucleotide polymorphism (SNP) ranking algorithm
 *
 *        Created:  08/11/2011
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu 
 *
 * =====================================================================================
 */

#ifndef __SNPRANK_H__
#define __SNPRANK_H__

#include <vector>
#include <string>
#include "armadillo"

using namespace std;
using namespace arma;

class SNPrank {

	public:
		SNPrank();
		SNPrank(string file);
		~SNPrank();
		void snprank(vector<string> names, mat  G, double gamma, string outFile);
		vector<string> getHeader();
		void setHeader(vector<string> hdr);
		mat getData();
		void setData(mat dat);
	private:
		void readFile(string filename);
		vector<string> header;
		mat data;
};
#endif

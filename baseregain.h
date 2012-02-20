/*
 * =====================================================================================
 *
 *       Filename:  baseregain.h
 *
 *    Description:  Base Regression GAIN class
 *
 *        Created:  02/19/2012
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu
 *
 * =====================================================================================
 */

#ifndef __BASEREGAIN_H__
#define __BASEREGAIN_H__

using namespace std;

class BaseRegain {
	public:
		virtual ~BaseRegain() {};
		virtual void run() = 0;
		virtual void writeRegain(bool fdrpr = false) = 0;
		virtual void writePvals() = 0;
		virtual void fdrPrune(double fdr) = 0;

};
#endif

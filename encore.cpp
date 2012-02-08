/*
 * =====================================================================================
 *
 *       Filename:  encore.cpp
 *
 *    Description:  Encore - A computational framework for analysis of GWAS and other 
 *                  biological data
 *
 *					Includes epistasis interaction analysis (reGAIN), eigenvector 
 *					centrality SNP ranking (SNPrank), and feature selection using the 
 *					Evaporative Cooling machine learning tool.  Encore also utilizes 
 *					plink for its ubiquitous data formats and wide array of GWAS 
 *					functionality.
 *
 *        Created:  02/01/2012 
 *
 *         Author:  Nick Davis, nick-davis@utulsa.edu
 *
 * =====================================================================================
 */

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "ec/EvaporativeCooling.h"
#include "plink/plinklibhandler.h"
#include "plink/plink.h"
#include "plink/options.h"
#include "plink/helper.h"

#include "snprank.h"
#include "regain.h"

using namespace boost;
namespace po = boost::program_options;

/********************************
 * required plink data structures
 *******************************/
// plink object
Plink* PP;

int main(int argc, char* argv[]) {
	// command line variables
	// data files
	string infile = "";
	string outfile_pref = "encore";
	string covarfile = "";
	string phenofile = "";
	string extrfile = "";
	//snprank
	double gamma = 0.85;	
	// reGAIN
	double sif_thresh = 0.05;
	double fdr = 0.5;
	// EC
	string ec_algo = "all";
	string ec_sm = "gm";
		   
	po::options_description desc("Encore - a tool for analysis of GWAS and other "
			"biological data.\nUsage:  encore -i snpdata.ped [mode] -o output-prefix");
	desc.add_options()
		("input-file,i", po::value<string>(&infile),
		 "Input GWAS file (.bed or .ped) or GAIN/reGAIN matrix (tab- or comma-separated)"
		)
		("output-prefix,o", po::value<string>(&outfile_pref),
		 "Prefix to use for all output files"
		)
		("snprank,s",
		 "Perform SNPrank analysis *mode*"
		)
			("gamma,g", po::value<double>(&gamma)->default_value(0.85, "0.85"),
			 "Damping factor"
			)
		("regain,r", 
		 "Calculate regression GAIN *mode*"
		)
			("compress-matrices", 
			 "Write binary (compressed) reGAIN matrices"
			)
			("sif-threshold", po::value<double>(&sif_thresh)->default_value(0.05, "0.05"),
			 "Numerical cutoff for SIF file interaction scores"
			)
			("fdr-prune", 
			 "FDR prune reGAIN interaction terms"
			)
			("fdr", po::value<double>(&fdr)->default_value(0.5, "0.5"),
			 "FDR value for BH method"
			)
		("ec,e", 
		 "Perform Evaporative Cooling (EC) analysis *mode*"
		)
			("ec-algorithm", po::value<string>(&ec_algo),
			 "EC ML algorithm (all|rf|rj)"
			)
			("ec-snp-metric", po::value<string>(&ec_sm),
			 "EC SNP metric (gm|am)"
			)
		("extract", po::value<string>(&extrfile),
		 "Extract list of SNPs from specified file"
		)
		("covar", po::value<string>(&covarfile),
		 "Include covariate file in analysis"
		)
		("pheno", po::value<string>(&phenofile),
		 "Include alternate phenotype file in analysis"
		)
		("assoc", 
		 "Run Case/control, QT association tests *mode*"
		)
		("linear", 
		 "Run linear regression model *mode*"
		)
		("ld-prune,l",
		 "Linkage disequilibrium (LD) pruning *mode*"
		)
		("help,h", 
		 "display this help screen"
		)
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    
	
	/********************************
	 * Help
	 *******************************/
	 if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	/*********************************
	 * Validate mode 
	 ********************************/
	if (!(vm.count("snprank") || vm.count("regain") || vm.count("ec") ||
				vm.count("assoc") || vm.count("linear") || vm.count("ld-prune"))) {
		cerr << "Error: Invalid command mode, must be one of:" << endl
			<< " --snprank, --regain, --ec, --assoc, --linear, --ld-prune" << endl
			<< endl << desc << endl;

		return 1;
	}

	/*********************************
	 * Validate only one mode passed
	 ********************************/
	int modes = 0;
	for (po::variables_map::iterator iter = vm.begin(); iter != vm.end(); ++iter) {
		if (iter->first == "snprank" || iter->first == "regain" ||
				iter->first == "ec" || iter->first == "assoc" ||
				iter->first == "linear" || iter->first == "ld-prune") {
					modes++;
				}
	}

	if (modes > 1) {
			cerr << "Error: Only one mode may be specified" << endl << endl << desc << endl;
			return 1;
	}

	/********************************
	 * Input file
	 *******************************/
	// require input file
	if (!vm.count("input-file")) {
		cerr << "Error: Must specify input file" << endl << endl << desc << endl;
		return 1;
	}

	// ensure SNP/reGAIN file exists
	else if (!boost::filesystem::exists(infile)) {
		cerr << "Error: Input file " << infile << " does not exist" << endl;
		return 1;
	}

	// Plink input file
	else if (infile.find(".bed") != string::npos || 
			infile.find(".ped") != string::npos) {
		// set file root
		vector<string> fileparts;
		split(fileparts, infile, is_any_of("."));
		par::fileroot = "";
		
		// handle files with multiple .s in the name
		for (int i =0; i < fileparts.size() - 1; i++) {
			par::fileroot += fileparts[i];
			if (fileparts.size() > 2 && i < fileparts.size() - 2) 
				par::fileroot += ".";
		}

		// set Plink's output file prefix
		par::output_file_name = outfile_pref;
		// initialize requisite Plink data structures
		initPlink();

		// read SNP file in PLINK
		// binary file
		if (infile.find(".bed") != string::npos) readPlBinFile(); 
		// plaintext file
		else if (infile.find(".ped") != string::npos) readPlFile();

		// additional PLINK setup
		initPlStats();
	}

	/********************************
	 * Covar file
	 *******************************/
	if (vm.count("covar")){
		// validate that covar file is used with proper modes
		if (!(vm.count("regain") || vm.count("linear"))) {
			cerr << "Error: Covariate file may only be used with --regain or --linear"
				<< endl << desc << endl;
			return 1;
		}

		// ensure covariate file exists
		else if (!boost::filesystem::exists(covarfile)) {
			cerr << "Error: Covariate file " << covarfile << " does not exist" << endl;
			return 1;
		}

		// read covariate file using PLINK
		else {
			par::covar_file = true;
			par::clist = true;
			par::clist_filename = covarfile;
			if(!PP->readCovListFile()){
				cerr << "Error: Problem reading the covariates" << endl;
				return 1;
			}
		}
	}
	
	/********************************
	 * Pheno file
	 *******************************/
	if (vm.count("pheno")) {
		// alternate phenotype validation
		if (vm.count("snprank")) {
			cerr << "Error: Alternate phenotype file cannot be used with "\
				"--snprank" << endl << desc << endl;
			return 1;
		}

		// ensure alternate phenotype file exists
		else if (!boost::filesystem::exists(phenofile)) {
			cerr << "Error: Alernate phenotype file " << phenofile << " does not exist" << endl;
			return 1;
		}

		// read alternate phenotype file using PLINK
		else {
			par::pheno_file = true;
			par::pheno_filename = phenofile;
			if(!PP->readPhenoFile())
				cerr << "Error: Problem reading the alternate phenotype file" << endl;
			}
	}

	/********************************
	 * Extract file
	 *******************************/
	if (vm.count("extract")) {
		// extract validation
		if (vm.count("snprank") || vm.count("ec") || vm.count("ldprune")) {
			cerr << "Error: Extract file cannot be used with "\
				"--snprank, --ec, or --ldprune" << endl << desc << endl;
			return 1;
		}

		// ensure extract file exists
		else if (!boost::filesystem::exists(extrfile)) {
			cerr << "Error: Extract file " << extrfile << " does not exist" << endl;
			return 1;
		}

		// read extract file SNPs using PLINK
		else {
			par::extract_set = true;
			par::extract_file = extrfile;
			PP->extractExcludeSet(false);
		}
	}
	
	/*********************************
	 * Validate mode sub-options
	 ********************************/
	if (!vm["gamma"].defaulted() && !vm.count("snprank")) {
			cerr << "Error: --gamma must be used with --snprank" << endl << endl <<
				desc << endl;
			return 1;
	}


	if (!vm["sif-threshold"].defaulted() && !vm.count("regain")) {
			cerr << "Error: --sif-threshold must be used with --regain" << endl << endl <<
				desc << endl;
			return 1;
	}

	if (vm.count("fdr-prune") && !vm.count("regain")) {
			cerr << "Error: --fdr-prune must be used with --regain" << endl << endl <<
				desc << endl;
			return 1;
	}

	if (!vm["fdr"].defaulted() && !vm.count("regain")) {
			cerr << "Error: --fdr must be used with --regain" << endl << endl <<
				desc << endl;
			return 1;
	}

	if (vm.count("ec-algorithm") && !vm.count("ec")) {
			cerr << "Error: ec-algorithm must be used with --ec" << endl << endl <<
				desc << endl;
			return 1;
	}

	if (vm.count("ec-snp-metric") && !vm.count("ec")) {
			cerr << "Error: ec-snp-metric must be used with --ec" << endl << endl <<
				desc << endl;
			return 1;
	}


	/*********************************
	 * Check primary mode of operation
	 ********************************/
	// SNPrank
	if (vm.count("snprank")) {
		SNPrank* sr = new SNPrank(infile);
		cout << "Writing SNPrank results to [ " << outfile_pref << ".snprank ]" << endl;
		sr->snprank(sr->getHeader(), sr->getData(), gamma, outfile_pref + ".snprank");
		delete sr;
	}

	// reGAIN
	else if (vm.count("regain")) {
		if (!vm.count("extract"))
			cout << "Warning: It is recommended to use an --extract file of SNPs with "\
				"--regain" << endl;
		// SNP major mode or individual major mode?
		if(par::fast_epistasis) {
			if(!par::SNP_major)
				PP->Ind2SNP();
		} else {
			if(par::SNP_major)
				PP->SNP2Ind();
		}

		bool fdrprune = vm.count("fdr-prune");
		Regain* r = new Regain(vm.count("compress-matrices"), sif_thresh, fdrprune);
		r->run();
		if (fdrprune){
			r->writeRegain();
			r->fdrPrune(fdr);
		}
		r->writeRegain(fdrprune);
		r->writePvals();
		delete r;
	}

	// Evaporative Cooling (EC)
	else if (vm.count("ec")) {
		// EC options map
		map<string,string> opts;
		// required options for EC
		opts.insert(pair<string,string>("ec-num-target", "0"));
		opts.insert(pair<string,string>("snp-data", infile));
		opts.insert(pair<string,string>("out-files-prefix", outfile_pref));

		// defaults for ID matching
		string numericfile = "";
		vector<string> ind_ids;
		vector<string> numeric_ids;
		vector<string> pheno_ids;
		bool datasetLoaded = false;

		// validate algorithm
		if (ec_algo != "all" && ec_algo != "rj" && ec_algo != "rf")
			cerr << "Error: EC algorithm must be one of: (all|rj|rf)" << endl;
		else opts.insert(pair<string,string>("ec-algorithm-steps", ec_algo));

		// validate metric
		if (ec_sm != "gm" && ec_sm != "am")
			cerr << "Error: EC SNP metric must be one of: (gm|am)" << endl;
		else opts.insert(pair<string,string>("snp-metric", ec_sm));

		// find IDs for loading from the dataset
		if(!GetMatchingIds(numericfile, phenofile,
						numeric_ids, pheno_ids, ind_ids))
			cerr << "Error: could not get matching IDs from numeric " <<
				"and/or phenotype files" << endl;
		// initialize dataset by extension
		Dataset* ds = 0;
		ds  = ChooseSnpsDatasetByExtension(infile);
		bool loaded = ds->LoadDataset(infile, "", phenofile, ind_ids);
		if (!loaded)
			cerr << "Error: Failure to load dataset for analysis" << endl;

		// file data stats 
		ds->PrintStats();

		// create ec object and run
		EvaporativeCooling* ec = new EvaporativeCooling(ds, opts, SNP_ONLY_ANALYSIS);
		if(!ec->ComputeECScores())
			cerr << "Error: Failed to calculate EC scores" << endl;

		// write results to file
		cout << "Writing EC results to [ " << outfile_pref << ".ec ]" << endl;
		ec->WriteAttributeScores(outfile_pref);
		delete ds;
		delete ec;
	}

	// Case/Control, QT association test OR linear model
	else if (vm.count("assoc")) {
		if (vm.count("linear"))
			par::assoc_glm = true;
		par::assoc_test = true;	
		PP->calcAssociationWithPermutation(*PP->pperm);
	}

	// LD-based pruning
	else if (vm.count("ldprune")) {
		par::prune_ld = true;
		par::prune_ld_pairwise = true;
		par::prune_ld_win = 50;
		par::prune_ld_step = 5;
		par::prune_ld_vif = 0.5;

		PP->pruneLD();
	}


	// Plink exit
	shutdown();
	return 0;
}



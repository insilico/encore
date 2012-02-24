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
#include <vector>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "ec/EvaporativeCooling.h"
#include "plink/plinklibhandler.h"
#include "plink/options.h"
#include "plink/helper.h"

#include "snprank.h"
#include "regain.h"

using namespace boost;
namespace po = boost::program_options;

/********************************
 * required plink data structures
 *******************************/

int main(int argc, char* argv[]) {
	// command line variables
	// data files
	string infile = "";
	string outfile_pref = "encore";
	string numfile = "";
	string covarfile = "";
	string phenofile = "";
	string extrfile = "";
	string remfile = "";
	string keepfile = "";
	// plink option defaults
	double maf = 0.0;
	double geno = 1;
	double mind = 1;
	double hwe = 0.001;
	double hwe2 = 0.001;
	//snprank
	double gamma = 0.85;	
	// reGAIN
	double sif_thresh = 0.05;
	double fdr = 0.5;
	// EC
	string ec_algo = "all";
	string ec_sm = "gm";
	double ec_numt = 0;
		   
	po::options_description desc("Encore - a tool for analysis of GWAS and other "
			"biological data.\nUsage:  encore -i snpdata.ped [mode] -o output-prefix");
	desc.add_options()
		("input-file,i", po::value<string>(&infile),
		 "Input GWAS file (.bed or .ped) or GAIN/reGAIN matrix (tab- or comma-separated)"
		)
		("output-prefix,o", po::value<string>(&outfile_pref),
		 "Prefix to use for all output files"
		)
		("numeric,n", po::value<string>(&numfile),
		 "Numeric file for quantitative data (uses PLINK covariate file format)"
		)
		("snprank,s",
		 "Perform SNPrank analysis *mode*"
		)
			("gamma", po::value<double>(&gamma)->default_value(0.85, "0.85"),
			 "SNPrank algorithm damping factor"
			)
		("regain,r", 
		 "Calculate regression GAIN *mode*"
		)
			("compress-matrices", 
			 "Write binary (compressed) reGAIN matrices"
			)
			("sif-threshold", po::value<double>(&sif_thresh)->default_value(0.05, "0.05"),
			 "Numerical cutoff for SIF file (generated by reGAIN) interaction scores"
			)
			("fdr-prune", 
			 "FDR prune reGAIN interaction terms"
			)
			("fdr", po::value<double>(&fdr)->default_value(0.5, "0.5"),
			 "FDR value for BH method applied to reGAIN"
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
			("ec-num-target", po::value<double>(&ec_numt),
			 "EC target number of attributes to keep"
			)
		("extract", po::value<string>(&extrfile),
		 "Extract list of SNPs from specified file"
		)
		("remove", po::value<string>(&remfile),
		 "Remove list of individuals from specified file"
		)
		("keep", po::value<string>(&keepfile),
		 "Keep list of individuals from specified file"
		)
		("prune",
		 "Remove individuals with missing phenotypes"
		)
		("covar", po::value<string>(&covarfile),
		 "Include covariate file in analysis"
		)
		("pheno", po::value<string>(&phenofile),
		 "Include alternate phenotype file in analysis"
		)
		("assoc", 
		 "Case/control, QTL association *mode*"
		)
		("linear", 
		 "Test for quantitative traits and multiple covariates *mode*"
		)
		("logistic",
		 "Test for disease traits and multiple covariates *mode*"
		)
		("model",
		 "Cochran-Armitage and full-model C/C association *mode*"
		)
		("model-trend",
		 "Use CA-trend test from model *mode*"
		)
		("model-gen",
		 "Use genotypic test from model *mode*"
		)
		("model-dom",
		 "Use dominant test from model *mode*"
		)
		("model-rec",
		 "Use recessive test from model *mode*"
		)
		("freq",
		 "Allele frequencies *mode*"
		)
		("counts",
		 "Modifies --freq to report actual allele counts"
		)
		("missing",
		 "Missing rates (per individual, per SNP) *mode*"
		)
		("maf", po::value<double>(&maf)->default_value(0.0, "0.0"),
		 "Minor allele frequency"
		)
		("geno", po::value<double>(&geno)->default_value(1, "1"),
		 "Maximum per-SNP missing"
		)
		("mind", po::value<double>(&mind)->default_value(1, "1"),
		 "Maximum per-person missing"
		)
		("hwe", po::value<double>(&hwe)->default_value(0.001, "0.001"),
		 "Hardy-Weinberg disequilibrium p-value (exact)"
		)
		("hwe2", po::value<double>(&hwe2)->default_value(0.001, "0.001"),
		 "Hardy-Weinberg disequilibrium p-value (asymptotic)"
		)
		("map3",
		 "Specify 3-column MAP file format"
		)
		("no-sex",
		 "PED file does not contain column 5 (sex)"
		)
		("allow-no-sex",
		 "Do not set ambiguously-sexed individuals missing"
		)
		("no-parents",
		 "PED file does not contain columns 3,4 (parents)"
		)
		("no-fid",
		 "PED file does not contain columns 1 (family ID)"
		)
		("r",
		 "Pairwise SNPxSNP LD (r)"
		)
		("r2",
		 "Pairwise SNPxSNP LD (r^2)"
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

	const char* margs[] = {
		"snprank",
		"regain",
		"ec",
		"assoc",
		"linear",
		"logistic",
		"model",
		"model-trend",
		"model-gen",
		"model-dom",
		"model-rec",
		"ld-prune",
		"freq",
		"missing",
		"r",
		"r2"
	};
	vector<string> modes(margs, margs + 16);
	PlinkHandler* ph;
	
	/********************************
	 * Help
	 *******************************/
	 if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	/*********************************
	 * Validate only one mode passed
	 ********************************/
	int nummodes = 0;
	for (po::variables_map::iterator iter = vm.begin(); iter != vm.end(); ++iter)
		if (find(modes.begin(), modes.end(), iter->first) != modes.end()) nummodes++;

	if (nummodes == 0) {
		cerr << "Error: Invalid command mode, must be one of:" << endl; 
		for (int i = 0; i < modes.size(); i++) 
			cerr << "\t--" << modes[i] << endl;
		cerr << endl << desc << endl;
		return 1;
	}

	else if (nummodes > 1) {
			cerr << "Error: Only one mode may be specified" << endl << endl << desc << endl;
			return 1;
	}
	
	/* Plink data file options *********************************/
	if (vm.count("map3"))
	     par::map3 = true;

	if (vm.count("no-sex"))
	     par::ped_skip_sex = true;

	if (vm.count("allow-no-sex"))
	     par::ignore_missing_sex = true;

	if (vm.count("no-parents"))
	     par::ped_skip_parents = true;

	if (vm.count("no-fid"))
		par::ped_skip_fid = true;

	/* Plink data file options *********************************/

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
		ph = new PlinkHandler();

		// read SNP file in PLINK
		// binary file
		if (infile.find(".bed") != string::npos) ph->readBedFile();
		// plaintext file
		else if (infile.find(".ped") != string::npos) ph->readPedFile();

	}

	else if (infile.find(".gain") == string::npos &&
			infile.find(".regain") == string::npos) {
		cerr << "Error:  Input file must be a PLINK data file (.ped/.bed)"
			<< " or a GAIN/reGAIN matrix (.gain/.regain)" << endl;
		return 1;
	}

	/********************************
	 * Allele frequencies
	 *******************************/
	if (vm.count("freq")) {
		par::af_write = true;

		// display MAF counts instead of freqs?
		if(vm.count("counts")) par::af_count = true;
	}


	/* Plink filtering options *********************************/
	 
	/********************************
	 * Missing rates 
	 *******************************/
	if (vm.count("missing")) par::report_missing = true;

	// Note:  Plink resets the following three values from their defaults 
	// in options.cpp of 0.01, 0.1, and 0.1, respectively, to 0.0, 1, and 1, 
	// respectively.  So we must always set these to match Plink's behavior.
	
	/********************************
	 * Minor allele frequency
	 *******************************/
	par::min_af = maf;

	/********************************
	 * Maximum per-SNP missing
	 *******************************/
	par::MAX_GENO_MISSING = geno;

	/********************************
	 * Maximum per-SNP missing
	 *******************************/
	par::MAX_IND_MISSING = mind;

	/********************************
	 * Hardy-Weinberg (exact)
	 *******************************/
	if (!vm["hwe"].defaulted()) {
		par::HWD_test = par::HWD_report = true;
		par::HWD_limit = hwe;
	}

	/********************************
	 * Hardy-Weinberg (asymptotic)
	 *******************************/
	if (!vm["hwe2"].defaulted()) {
		par::HWD_test = par::HWD_report = true;
		par::HWD_standard = true;
		par::HWD_limit = hwe2;
	}

	/* end Plink filtering options ******************************/

	// additional PLINK setup
	if (!vm.count("snprank"))
		ph->initData();

	/********************************
	 * Numeric file
	 *******************************/
	// validate that numeric file is used with proper modes
	if (vm.count("numeric")) {
		if(!vm.count("regain")) {
			cerr << "Error: Numeric file may only be used with --regain"
				<< endl << desc << endl;
			return 1;
		}

		// ensure numeric file exists
		else if (!boost::filesystem::exists(numfile)) {
			cerr << "Error: Numeric file " << numfile << " does not exist" << endl;
			return 1;
		}

		// read the numeric attributes using PLINK
		else ph->readNumFile(numfile);

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
		else ph->readCovFile(covarfile);
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
		else ph->readPhenoFile(phenofile);
	}

	/********************************
	 * Prune
	 *******************************/
	if (vm.count("prune")) {
		// prune validation
		if (vm.count("snprank") || vm.count("ec") || vm.count("ld-prune")) {
			cerr << "Error: prune file cannot be used with "\
				"--snprank, --ec, or --ld-prune" << endl << desc << endl;
			return 1;
		}

		ph->pruneInd();
	}

	/********************************
	 * Extract file
	 *******************************/
	if (vm.count("extract")) {
		// extract validation
		if (vm.count("snprank") || vm.count("ec") || vm.count("ld-prune")) {
			cerr << "Error: Extract file cannot be used with "\
				"--snprank, --ec, or --ld-prune" << endl << desc << endl;
			return 1;
		}

		// ensure extract file exists
		else if (!boost::filesystem::exists(extrfile)) {
			cerr << "Error: Extract file " << extrfile << " does not exist" << endl;
			return 1;
		}

		// read extract file SNPs using PLINK
		else ph->readExtractFile(extrfile);
	}
	
	/********************************
	 * Remove file
	 *******************************/
	if (vm.count("remove")) {
		// remove validation
		if (vm.count("snprank") || vm.count("ec") || vm.count("ld-prune")) {
			cerr << "Error: Remove file cannot be used with "\
				"--snprank, --ec, or --ld-prune" << endl << desc << endl;
			return 1;
		}

		// ensure remove file exists
		else if (!boost::filesystem::exists(remfile)) {
			cerr << "Error: Remove file " << remfile << " does not exist" << endl;
			return 1;
		}

		// read remove file individuals using PLINK
		else ph->readRemoveFile(remfile);
	}

	/********************************
	 * Keep file
	 *******************************/
	if (vm.count("keep")) {
		// keep validation
		if (vm.count("snprank") || vm.count("ec") || vm.count("ld-prune")) {
			cerr << "Error: keep file cannot be used with "\
				"--snprank, --ec, or --ld-prune" << endl << desc << endl;
			return 1;
		}

		// ensure keep file exists
		else if (!boost::filesystem::exists(keepfile)) {
			cerr << "Error: keep file " << keepfile << " does not exist" << endl;
			return 1;
		}

		// read keep file individuals using PLINK
		else ph->readKeepFile(keepfile);
	}


	/*********************************
	 * Validate mode sub-options
	 ********************************/
	if (!vm["gamma"].defaulted() && !vm.count("snprank")) {
			cerr << "Error: --gamma must be used with --snprank" << endl << endl <<
				desc << endl;
			return 1;
	}

	if (vm.count("compress-matrices") && !vm.count("regain")) {
			cerr << "Error: --compress-matrices must be used with --regain" << endl << endl <<
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

	if (vm.count("ec-num-target") && !vm.count("ec")) {
			cerr << "Error: ec-num-target must be used with --ec" << endl << endl <<
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

		// set individual major mode
		ph->setInd();
		bool fdrprune = vm.count("fdr-prune");
		Regain* r;
		// if --numeric passed, create an integrative regain object
		if (vm.count("numeric")) r = new Regain(vm.count("compress-matrices"), sif_thresh, true,  fdrprune);
		else 
			r = new Regain(vm.count("compress-matrices"), sif_thresh, false, fdrprune);
		r->run();
		if (fdrprune){
			r->writeRegain(false);
			r->fdrPrune(fdr);
		}
		r->writeRegain(false, fdrprune);
		r->writeRegain(true);
		delete r;
	}

	// Evaporative Cooling (EC)
	else if (vm.count("ec")) {
		// EC options map
		map<string,string> opts;
		// required options for EC
		opts.insert(pair<string,string>("ec-num-target", lexical_cast<string>(ec_numt)));
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

	// Association tests and models
	else if (vm.count("assoc") || vm.count("linear") || vm.count("logistic") ||
			vm.count("model") || vm.count("trend") ||
			vm.count("model-trend") || vm.count("model-gen") ||
			vm.count("model-dom") || vm.count("model-rec")) {

		par::assoc_test = true;	

		if (vm.count("model") || vm.count("trend") ||
			vm.count("model-dom") || vm.count("model-rec") ||
			vm.count("model-trend") || vm.count("model-gen"))
				par::full_model_assoc = true;

		if (vm.count("trend")) {
			par::trend_only = true;
			par::model_perm_trend = true;
		}
		else {
			if (vm.count("model-gen")) par::model_perm_gen = true;
			else if (vm.count("model-dom")) par::model_perm_dom = true;
			else if (vm.count("model-rec")) par::model_perm_rec = true;
			else if (vm.count("model-trend")) par::model_perm_trend = true;
			else par::model_perm_best = true;
		}

		if (vm.count("linear") || vm.count("logistic")) par::assoc_glm = true;

		ph->assocTest();
	}

	// LD-based pruning
	else if (vm.count("ld-prune")) {
		ph->LDPrune();
	}

	// LD-pruning using r and r^2
  	else if(vm.count("r") || vm.count("r2")) {
  		if (vm.count("r2")) par::disp_r2 = true;
		else if (vm.count("r")) {
			par::disp_r1 = true;
			par::disp_r2 = false;
		}
		
		ph->LDStats();
	}

	// Plink exit
	if (!vm.count("snprank"))
		shutdown();
	return 0;
}



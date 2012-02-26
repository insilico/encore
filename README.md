Encore
========

#### A computational tool for integrative analysis of GWAS and other biological data ####

### Description ###
Encore is a free, open-source command-line tool for analysis of GWAS (SNP) and
other types of biological data.  Several modes are available for various types
of analysis, including:

 * Regression of SNPs and quantitative data in epistasis interaction networks
 (Regression [GAIN](http://insilico.utulsa.edu/gain)/reGAIN)
 * Eigenvector centrality ranking of top SNPs/quantitative attributes
 ([SNPrank](http://insilico.utulsa.edu/snprank))
 * Machine learning feature selection algorithms useful for filtering initially
   large data sets to the top few thousands attributes for subsequent analysis
 ([Evaporative Cooling](http://insilico.utulsa.edu/evaporative-cooling)/EC)
 * Data formats and encoding of SNP and quantitative data types
 (PLINK)
 * Several additional options ported from the PLINK library

Many libraries are used to provide various functionality.  The most
prominent is PLINK, a ubiquitous third-party GWAS analysis tool.  We have
modified the latest stable source release of PLINK into a library used to handle
GWAS formats, filter data, and provide statistical/association tests.  See our
PLINK project page for more details.

Evaporative Cooling (EC) is another key library, and provides feature selection
of SNPs and quantitative data, using ReliefF and Random Jungle for interactions
and main effects, respectively.  EC is also available as a standalone
[tool](http://insilico.utulsa.edu/evaporative-cooling).

Encore is being developed by the In Silico Research Group at the Tandy School
of Computer Science of the [University of Tulsa](http://www.utulsa.edu).  Our
research is sponsored by the NIH and William K. Warren foundation.  For more
details, visit our research [website](http://insilico.utulsa.edu).

### Dependencies ###
* EC library, available as a source release on the EC project page, and its
dependencies:

  * [Random Jungle](http://github.com/insilico/randomjungle)

  * gfortran, sometimes installed alongside compiler tools

  * GNU Scientific library (libgsl)

  * libxml2

* [Boost](http://www.boost.org) system, filesystem, and program-options libraries 
 
* The libz/zlib compression library is required, but this is installed by default
on most Unix systems.  In MinGW libz is installed via mingw-get.

* [Armadillo](http://arma.sourceforge.net), LAPACK and BLAS for the optimized
linear algebra code

* OpenMP is required to take advantage of the parallelized regression GAIN 
analysis code.  This is another library typically installed alongside the 
compiler toolchain.

### Compilation Environment and Instructions ###
To compile this code, a GNU toolchain and suitable environment are required.
GNU g++ has been used to successfully compile the code.

We have successfully built and run Encore on:

 * Linux (64-bit Ubuntu) (gcc-4.6)
 * Mac (10.6 - 10.7) (gcc-4.2.1)
 * Windows 7 (32-bit) using the [MinGW](http://www.mingw.org) compiler system
  (gcc-4.6)

To build Encore, first run the bootstrap script

    ./bootstrap.sh

This calls autoreconf and generates the configure script.  From this point, a 
standard

    ./configure && make && sudo make install

will generate the `Makefile`, compile and link the code, and copy the objects to
the installation directory (default of `/usr/local`).  As is convention, headers
are installed in `$PREFIX/include`, binary in `$PREFIX/bin`, and the library in
`$PREFIX/lib`.

The resulting binary encore/encore.exe will run as a command-line tool.

### Usage ###
To see the help screen,

    ./encore -h

    Encore - a tool for analysis of GWAS and other biological data.
    Usage:  encore -i snpdata.ped [mode] -o output-prefix:
      -i [ --input-file ] arg       Input GWAS file (.bed or .ped) or GAIN/reGAIN 
                                    matrix (tab- or comma-separated)
      -o [ --output-prefix ] arg    Prefix to use for all output files
      -n [ --numeric ] arg          Numeric file for quantitative data (uses PLINK 
                                    covariate file format)
      -d [ --data-summary ]         Simply print input file stats (for PLINK 
                                    .ped/.bed files)
      -s [ --snprank ]              Perform SNPrank analysis *mode*
      --gamma arg (=0.85)           SNPrank algorithm damping factor
      -r [ --regain ]               Calculate regression GAIN *mode*
      --compress-matrices           Write binary (compressed) reGAIN matrices
      --sif-threshold arg (=0.05)   Numerical cutoff for SIF file (generated by 
                                    reGAIN) interaction scores
      --fdr-prune                   FDR prune reGAIN interaction terms
      --fdr arg (=0.5)              FDR value for BH method applied to reGAIN
      -e [ --ec ]                   Perform Evaporative Cooling (EC) analysis 
                                    *mode*
      --ec-algorithm arg            EC ML algorithm (all|rf|rj)
      --ec-snp-metric arg           EC SNP metric (gm|am)
      --ec-num-target arg           EC target number of attributes to keep
      --extract arg                 Extract list of SNPs from specified file
      --remove arg                  Remove list of individuals from specified file
      --keep arg                    Keep list of individuals from specified file
      --prune                       Remove individuals with missing phenotypes
      --covar arg                   Include covariate file in analysis
      --pheno arg                   Include alternate phenotype file in analysis
      --make-bed                    Make .bed, .fam and .bim
      --ci arg (=0.95)              Confidence interval for CMH odds ratios
      --assoc                       Case/control, QTL association *mode*
      --linear                      Test for quantitative traits and multiple 
                                    covariates *mode*
      --logistic                    Test for disease traits and multiple covariates
                                    *mode*
      --model                       Cochran-Armitage and full-model C/C association
                                    *mode*
      --model-trend                 Use CA-trend test from model *mode*
      --model-gen                   Use genotypic test from model *mode*
      --model-dom                   Use dominant test from model *mode*
      --model-rec                   Use recessive test from model *mode*
      --freq                        Allele frequencies *mode*
      --counts                      Modifies --freq to report actual allele counts
      --missing                     Missing rates (per individual, per SNP) *mode*
      --missing-genotype arg (="0") Missing genotype code
      --maf arg (=0.0)              Minor allele frequency
      --geno arg (=1)               Maximum per-SNP missing
      --mind arg (=1)               Maximum per-person missing
      --hwe arg (=0.001)            Hardy-Weinberg disequilibrium p-value (exact)
      --hwe2 arg (=0.001)           Hardy-Weinberg disequilibrium p-value 
                                    (asymptotic)
      --filter-founders             Include only founders
      --map3                        Specify 3-column MAP file format
      --no-sex                      PED file does not contain column 5 (sex)
      --allow-no-sex                Do not set ambiguously-sexed individuals 
                                    missing
      --no-parents                  PED file does not contain columns 3,4 (parents)
      --no-fid                      PED file does not contain columns 1 (family ID)
      --r                           Pairwise SNPxSNP LD (r)
      --r2                          Pairwise SNPxSNP LD (r^2)
      -l [ --ld-prune ]             Linkage disequilibrium (LD) pruning *mode*
      -h [ --help ]                 display this help screen

All commands will include an input file (`-i/--input-file`), and, optionally, 
an output file prefix (`-o/--output-prefix`).  Each command must include a 
single mode, denoted in the help screen with `*mode*`.

To perform a regression GAIN (reGAIN) analysis,

    ./encore -i snpdata.ped --regain -o result

This will use genotype/phenotype information from `snpdata.ped`, a PLINK
plaintext GWAS file, in the regression.  All of the output files produced will
be prepended with 'result'.

Once a `.regain` file is produced, this can be used as input for eigenvector 
centrality ranking with SNPrank.

    ./encore -i result.regain --snprank -o rankings


This produces a file called `rankings.snprank`, in which the SNPs are ranked 
in descending order, along with their associated SNPrank and information gain
(IG) scores.

To invoke EC to filter the original SNP data,

    ./encore -i snpdata.ped --ec -o filtered

Additional EC-based options can be specified, such as `--ec-algorithm` and
`--ec-snp-metric` to select the EC algorithm and SNP metric, respectively.

For additional examples, see the [Encore](http://insilico.utulsa.edu/encore)
page on our research website.

### Contributors ###
See AUTHORS file.

### References ###
N.A. Davis, J.E. Crowe, Jr., N.M. Pajewski, and B.A. McKinney. Surfing a
genetic association interaction network to identify modulators of antibody
response to smallpox vaccine. Genes and Immunity, 2010,
doi: 10.1038/gene.2010.3.

B.A. McKinney, J.Guo, J.E. Crowe, Jr., and D. Tian. Capturing the spectrum of
interaction effects in genetic association studies by simulated evaporative
cooling network analysis. PLoS Genetics 2009, 5(3): e1000432.
doi:10.1371/journal.pgen.1000432.

B.A. McKinney, D. M. Reif, B. C. White, J. E. Crowe Jr., J. H. Moore.
"Evaporative cooling feature selection for genotypic data involving
interactions." Bioinformatics 2007, 23: 2113-2120

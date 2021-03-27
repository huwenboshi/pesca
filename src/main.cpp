#include "utils.h"
#include "mcmc.h"
#include "numeric.h"
#include "fit.h"
#include "post.h"
#include <time.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <set>

using namespace std;
using namespace Eigen;

int main(int ac, char* av[]) {

    // set random seed
    srand(62442); 
    
    // parse command line
    po::variables_map vm = get_command_line(ac, av);
   
    string zsc_file_nm1, zsc_file_nm2, ld1_nm, ld2_nm, mode, out;
    double sigmasq1, sigmasq2, tot_nsnp;
    try {
        zsc_file_nm1 = vm["zscore1"].as<string>();
        zsc_file_nm2 = vm["zscore2"].as<string>();
        ld1_nm = vm["ld1"].as<string>();
        ld2_nm = vm["ld2"].as<string>();
        mode = vm["mode"].as<string>();
        out = vm["out"].as<string>();
        sigmasq1 = vm["sigmasq1"].as<double>();
        sigmasq2 = vm["sigmasq2"].as<double>();
        tot_nsnp = vm["totnsnp"].as<double>();
    }
    catch(...) {
        cerr <<"The following flags must be specified:" << endl
             <<"\t--mode\t\t(\"fit\" or \"post\")" << endl
             <<"\t--zscore1\t(a list of z-score files for pop 1)" << endl
             <<"\t--zscore2\t(a list of z-score files for pop 2)" << endl
             <<"\t--ld1\t\t(a list of LD files for pop 1)" << endl
             <<"\t--ld2\t\t(a list of LD files for pop 2)" << endl
             <<"\t--sigmasq1\t(heritability times sample size for pop 1)"<<endl
             <<"\t--sigmasq2\t(heritability times sample size for pop 2)"<<endl
             <<"\t--totnsnp\t(total number of SNPs)" << endl
             <<"\t--out\t\t(output file name)" << endl;
        cerr <<"The following flags are optional:" << endl
             <<"\t--nchain\t(number of MCMC chains, default is 1)" << endl
             <<"\t--templ\t\t(lower bound of simulated"
             <<" annealing temperature, default is 1.0)" << endl
             <<"\t--temph\t\t(higher bound of simulated"
             <<" annealing temperature, default is 1.0)" << endl
             <<"\t--nburn\t\t(number of MCMC burn-ins, default is 20000)"<<endl
             <<"\t--nsample\t(number of MCMC samples, default is 30000)"<<endl
             <<"\t--dist\t\t(\"mvb\" or \"mult\", default is \"mult\")"<<endl
             <<"\t--f00\t\t(f00 parameter of the MVB, default is 0.0)"<<endl
             <<"\t--f01\t\t(f01 parameter of the MVB, default is -6.89)"<<endl
             <<"\t--f10\t\t(f10 parameter of the MVB, default is -6.89)"<<endl
             <<"\t--f11\t\t(f11 parameter of the MVB, default is 8.98)"<<endl
             <<"\t--p00\t\t(p00 parameter of the Mult, default is 0.99)"<<endl
             <<"\t--p01\t\t(p01 parameter of the Mult, default is 0.001)"<<endl
             <<"\t--p10\t\t(p10 parameter of the Mult, default is 0.001)"<<endl
             <<"\t--p11\t\t(p11 parameter of the Mult, default is 0.008)"<<endl
             <<"\t--lambda\t(shrinkage parameter, default is 0.0001)"<<endl
             <<"\t--max_iter_fit\t(max number of EM iterations for prior"
             <<" estimation, default is 100)"<<endl
             <<"\t--max_iter_post\t(number of iterations for posterior"
             <<" estimation, default is 1)"<<endl
            <<"\t--print\t\t(\"yes\" or \"no\", default is \"yes\")" << endl;
        return 1;
    }
    
    double templ = vm["templ"].as<double>();
    double temph = vm["temph"].as<double>();
    size_t nchain = vm["nchain"].as<size_t>();
    size_t nburn = vm["nburn"].as<size_t>();
    size_t nsample = vm["nsample"].as<size_t>();
    size_t max_iter_fit = vm["max_iter_fit"].as<size_t>();
    size_t max_iter_post = vm["max_iter_post"].as<size_t>();
    double f00 = vm["f00"].as<double>();
    double f01 = vm["f01"].as<double>();
    double f10 = vm["f10"].as<double>();
    double f11 = vm["f11"].as<double>();
    double p00 = vm["p00"].as<double>();
    double p01 = vm["p01"].as<double>();
    double p10 = vm["p10"].as<double>();
    double p11 = vm["p11"].as<double>();
    double lambda = vm["lambda"].as<double>();
    string print = vm["print"].as<string>();
    string dist = vm["dist"].as<string>();

    // create files to write
    ofstream outfile, logfile;
    if(mode == "fit") {
        logfile.open(out+".log");
    }
    if(mode == "post") {
        outfile.open(out+".txt");
    }

    // load list
    vector<string> zscore1_list = load_list(zsc_file_nm1.c_str());
    vector<string> zscore2_list = load_list(zsc_file_nm2.c_str());
    vector<string> ld1_list = load_list(ld1_nm.c_str());
    vector<string> ld2_list = load_list(ld2_nm.c_str());
    size_t nloci = zscore1_list.size();

    // load list of zscores, ld, and causal effect size
    vector<vector<VectorXf> > all_zsc;
    vector<vector<MatrixXf> > all_ld_mat;
    vector<vector<VectorXf> > all_caueff;
    vector<vector<VectorXf> > all_caueff_var;
    vector<vector<double> > all_sigma_sq;
    vector<vector<Zscore_t> > all_zsc_inf1;
    vector<vector<Zscore_t> > all_zsc_inf2;
    vector<size_t> all_nsnp;
    for(size_t i=0; i<nloci; i++) {
  
        // load z scores and ld matrices
        vector<Zscore_t> zscore_info1 = load_zscore(zscore1_list[i].c_str());
        vector<Zscore_t> zscore_info2 = load_zscore(zscore2_list[i].c_str());
        size_t nsnp = zscore_info1.size();
        VectorXf zsc1 = get_zsc_vec(zscore_info1);
        VectorXf zsc2 = get_zsc_vec(zscore_info2);
        
        MatrixXf ld1 = load_ldmat(ld1_list[i].c_str());
        MatrixXf ld2 = load_ldmat(ld2_list[i].c_str());

        // record number of snps
        all_nsnp.push_back(nsnp);

        // create a vector
        all_zsc.push_back(vector<VectorXf>());
        all_ld_mat.push_back(vector<MatrixXf>());

        // append the data
        all_zsc[i].push_back(zsc1);
        all_ld_mat[i].push_back(ld1);
        all_zsc_inf1.push_back(zscore_info1);
        all_zsc[i].push_back(zsc2);
        all_ld_mat[i].push_back(ld2);
        all_zsc_inf2.push_back(zscore_info2);
    }

    // find the sigma to use
    for(size_t i=0; i<all_nsnp.size(); i++) {
        double frac = ((double)all_nsnp[i])/tot_nsnp;
        double sigmasq1_locus = frac*sigmasq1;
        double sigmasq2_locus = frac*sigmasq2;
        all_sigma_sq.push_back(vector<double>());
        all_sigma_sq[i].push_back(sigmasq1_locus);
        all_sigma_sq[i].push_back(sigmasq2_locus);
    }

    // initialize parameter
    VectorXf params = VectorXf::Zero(4);
    if(dist == "mvb") {
        params(0) = f00;
        params(1) = f01;
        params(2) = f10;
        params(3) = f11;
    }
    else if(dist == "mult") {
        if(p01 <= 0.0 || p01 <= 0.0 || p10 <= 0.0 || p11 <= 0.0) {
            cerr << "Parameters cannot be negative for Mult" << endl;
            return 1;
        }
        params(0) = 0.0;
        params(1) = log(p01/p00);
        params(2) = log(p10/p00);
        params(3) = log(p11/p00) - params(2) - params(1);
    }
    else {
        cerr << "Unrecognized distribution" << endl;
        return 1;
    }

    // perform fitting
    if(mode == "fit" && !all_nsnp.empty()) {
        fit(all_zsc, all_sigma_sq, all_ld_mat, all_nsnp, params, templ,
            temph, nchain, nburn, nsample, lambda, max_iter_fit,
            outfile, logfile, print);
    }

    // perform posterior inference
    else if(mode == "post" && !all_nsnp.empty()) {
        post(all_zsc, all_sigma_sq, all_ld_mat, all_nsnp, params, templ,
            temph, nchain, nburn, nsample, lambda, max_iter_post, all_zsc_inf1,
            all_zsc_inf2, outfile, logfile);
    }

    // close files
    if(mode == "post") {
        outfile.close();
    }
    if(mode == "fit") {
        logfile.close();
    }
}

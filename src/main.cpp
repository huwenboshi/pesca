#include "mcmc.h"
#include "numeric.h"
#include "utils.h"
#include "fit.h"
#include "post.h"
#include <time.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <set>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(int ac, char* av[]) {

    // set random seed
    srand(time(NULL)); 
    
    // parse command line
    po::variables_map vm = get_command_line(ac, av);
    string zsc_file_nm1 = vm["zscore1"].as<string>();
    string zsc_file_nm2 = vm["zscore2"].as<string>();
    string ld1_nm = vm["ld1"].as<string>();
    string ld2_nm = vm["ld2"].as<string>();
    string mode = vm["mode"].as<string>();
    string out = vm["out"].as<string>();
    double sigmasq1 = vm["sigmasq1"].as<double>();
    double sigmasq2 = vm["sigmasq2"].as<double>();
    double tot_nsnp = vm["totnsnp"].as<double>();
    double templ = vm["templ"].as<double>();
    double temph = vm["temph"].as<double>();
    size_t nchain = vm["nchain"].as<size_t>();
    size_t nburn = vm["nburn"].as<size_t>();
    size_t nsample = vm["nsample"].as<size_t>();
    size_t max_iter = vm["max_iter"].as<size_t>();
    double f00 = vm["f00"].as<double>();
    double f01 = vm["f01"].as<double>();
    double f10 = vm["f10"].as<double>();
    double f11 = vm["f11"].as<double>();
    double lambda = vm["lambda"].as<double>();

    // create files to write
    ofstream outfile, logfile;
    outfile.open(out+".txt");
    logfile.open(out+".log");

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

    // perform fitting
    if(mode == "fit" && !all_nsnp.empty()) {
        
        // initialize parameter
        VectorXf params = VectorXf::Zero(4);
        params(0) = f00;
        params(1) = f01;
        params(2) = f10;
        params(3) = f11;

        // estimate the parameters
        fit(all_zsc, all_sigma_sq, all_ld_mat, all_nsnp, params, templ,
            temph, nchain, nburn, nsample, lambda, max_iter, outfile, logfile);
    }

    // perform posterior inference
    else if(mode == "post" && !all_nsnp.empty()) {
       
        // initialize parameter
        VectorXf params(4);
        params(0) = f00;
        params(1) = f01;
        params(2) = f10;
        params(3) = f11;

        // estimate the parameters
        post(all_zsc, all_sigma_sq, all_ld_mat, all_nsnp, params, templ,
            temph, nchain, nburn, nsample, lambda, max_iter, all_zsc_inf1,
            all_zsc_inf2, outfile, logfile);
    }

    // close files
    outfile.close();
    logfile.close();
}

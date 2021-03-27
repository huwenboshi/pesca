#include "fit.h"
#include "numeric.h"

using namespace std;

#define DEBUG_FIT true

// fit the model to obtain prior estimation
void fit(const vector<vector<VectorXf> >& all_zsc,
         const vector<vector<double> >& all_sigma_sq,
         const vector<vector<MatrixXf> >& all_ld,
         const vector<size_t>& all_nsnp, VectorXf& params,
         double templ, double temph, size_t nchain,
         size_t nburn, size_t nsample, double lambda,
         size_t max_iter, ofstream& outfile, ofstream& logfile,
         string print) {

    // count total number of SNPs in the loci
    double nsnp_tot = 0.0;
    for(size_t l=0; l<all_nsnp.size(); l++) {
        nsnp_tot += all_nsnp[l];
    }
    size_t nloci = all_zsc.size();
   
    // write initialization
    logfile << "iter\tnsnp\tq00\tq01\tq10\tq11\tf00\tf01\tf10\tf11" << endl;
    double denom = log_mvb_denom(params);
    logfile << "0" << "\t" << nsnp_tot << "\t";
    if(print == "yes") {
        cerr << "iter\tnsnp\tq00\tq01\tq10\tq11\tf00\tf01\tf10\tf11" << endl;
        cerr << "0" << "\t" << nsnp_tot << "\t";
    }
    
    for(size_t i=0; i < params.rows(); i++) {
        logfile << nsnp_tot*exp(mvb_sfunc(params, i)-denom) << "\t";
        if(print == "yes") {
            cerr << nsnp_tot*exp(mvb_sfunc(params, i)-denom) << "\t";
        }
    }
    logfile << "\t";
    if(print == "yes") {
        cerr << "\t";
    }
    for(size_t i=0; i < params.rows(); i++) {
        logfile << params(i);
        if(print == "yes") {
            cerr << params(i);
        }
        if(i < params.rows() - 1) {
            logfile << "\t";
            if(print == "yes") {
                cerr << "\t";
            }
        }
    }
    logfile << endl;
    if(print == "yes") {
        cerr << endl;
    }

    // create one sampler for each locus
    double tstep = (temph-templ)/(((double)nchain)-1+eps);
    vector<Sampler> all_sampler;
    for(size_t l=0; l<nloci; l++) {
        std::vector<Chain> chains;
        for(size_t i=0; i<nchain; i++) {
            double temp = templ + i*tstep;
            chains.push_back(Chain(all_zsc[l], all_sigma_sq[l],
                all_ld[l], params, temp, lambda));
        }
        Sampler sampler(chains, nburn, nsample);
        all_sampler.push_back(sampler);
    }

    // iterate until converge
    for(size_t iter=0; iter<max_iter; iter++) {

        // obtain the count
        VectorXf conf_cnt = VectorXf::Zero(params.rows());
        for(size_t l=0; l<nloci; l++) {
            VectorXf conf_cnt_locus = VectorXf::Zero(params.rows());
            while(!all_sampler[l].done()) {

                // get a sample
                vector<config_t> sample = all_sampler[l].sample();
                config_t conf1 = sample[0];
                config_t conf2 = sample[1];

                // count the occurance of configurations
                double n11 = 0.0;
                for(auto it=conf1.begin(); it != conf1.end(); it++) {
                    n11 += (conf2.find(*it) != conf2.end());
                }
            
                // record the count
                conf_cnt_locus(0) += all_nsnp[l];
                conf_cnt_locus(1) += conf1.size();
                conf_cnt_locus(2) += conf2.size();
                conf_cnt_locus(3) += n11;
            }
            for(size_t i=0; i < conf_cnt_locus.rows(); i++) {
                conf_cnt_locus(i) = conf_cnt_locus(i)/nsample;
            }
            conf_cnt += conf_cnt_locus;
        }
        
        // using gradient ascent to complete the m step
        update_params(params, conf_cnt);

        // reset the sampler with new parameters
        for(size_t l=0; l < nloci; l++) {
            all_sampler[l].reset(params);
        }

        // log the output of each iteration
        double denom = log_mvb_denom(params);
        logfile << iter+1 << "\t" << nsnp_tot << "\t";
        if(print == "yes") {
            cerr << iter+1 << "\t" << nsnp_tot << "\t";
        }
        for(size_t i=0; i < params.rows(); i++) {
            logfile << nsnp_tot*exp(mvb_sfunc(params, i)-denom) << "\t";
            if(print == "yes") {
                cerr << nsnp_tot*exp(mvb_sfunc(params, i)-denom) << "\t";
            }
        }
        logfile << "\t";
        if(print == "yes") {
            cerr << "\t";
        }
        for(size_t i=0; i < params.rows(); i++) {
            logfile << params(i);
            if(print == "yes") {
                cerr << params(i);
            }
            if(i < params.rows() - 1) {
                logfile << "\t";
                if(print == "yes") {
                    cerr << "\t";
                }
            }
        }
        logfile << endl;
        if(print == "yes") {
            cerr << endl;
        }
    }
}

// update parameters
void update_params(VectorXf& params, const VectorXf& conf_cnt) {

    double pdcnt = 0.001;

    double n00 = conf_cnt(0)-conf_cnt(1)-conf_cnt(2)+conf_cnt(3) + pdcnt;
    double n01 = conf_cnt(1)-conf_cnt(3) + pdcnt;
    double n10 = conf_cnt(2)-conf_cnt(3) + pdcnt;
    double n11 = conf_cnt(3) + pdcnt;

    params(0) = 0.0;
    params(1) = log(n01)-log(n00);
    params(2) = log(n10)-log(n00);
    params(3) = log(n11)-log(n00)-params(1)-params(2);
}

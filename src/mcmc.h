#include <set>
#include <string>
#include "utils.h"
#include "numeric.h"
#include <stdlib.h>
#include <numeric>
#include <algorithm>

using namespace std;

#ifndef MCMC
#define MCMC

// the mcmc chain
class Chain {

    public:
        Chain(const vector<VectorXf>& zsc,
              const vector<double>& sigma_sq,
              const vector<MatrixXf>& ld,
              const VectorXf& params,
              double temp, double lamba);
        void next_state();
        vector<config_t> get_state();
        void set_state(const vector<config_t>& state);
        void reset(const VectorXf& params);
        vector<VectorXf> get_zsc();
        double get_loglik_mvb();
        double eval_loglik_mvb(vector<config_t> config);
        double get_temp();
        
    private:
        const vector<VectorXf>& m_zsc;
        const vector<double>& m_sigma_sq;
        const vector<MatrixXf>& m_ld;
        VectorXf m_params;
        vector<config_t> m_config;
        vector<double> m_logbf;
        double log_bf(size_t updt, const config_t& conf);
        double m_loglik_mvb;
        double m_temp;
        double m_lambda;
        size_t m_nsnp;
        size_t m_nstep;
        size_t m_npop;
};

// the sampler that generate samples
class Sampler {

    public:
        Sampler(const vector<Chain>& chains, size_t nburn, size_t nsample);
        void reset(const VectorXf& params);
        vector<config_t> sample();
        bool done();

    private:
        vector<Chain> m_chains;
        size_t m_nburn;
        size_t m_nsample;
        size_t m_ndone;
        void swap_chain();
};

// generate a random integer between min and max inclusive
size_t rand_int(size_t min, size_t max);

// generate a random number between 0 and 1
double rand_double();

#endif

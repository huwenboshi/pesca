#include "mcmc.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <random>
#include <cstdlib>

using namespace std;

// generate a random integer between min and max inclusive
size_t rand_int(size_t min, size_t max) {
    return rand()%(max-min+1)+min;
}

// generate a random number between 0 and 1
double rand_double() {
    return rand() / (RAND_MAX + 1.0);
}

// initialize the chain
Chain::Chain(const vector<VectorXf>& zsc,
             const vector<double>& sigma_sq,
             const vector<MatrixXf>& ld,
             const VectorXf& params,
             double temp, double lambda):
        m_zsc(zsc), m_sigma_sq(sigma_sq), m_ld(ld),
        m_params(params), m_temp(temp), m_lambda(lambda){
    
    // initialize the counter
    m_nsnp = ld[0].rows();
    m_nstep = 0;

    // initialize the configs, mean, and loglik for each zscore
    m_npop = zsc.size();
    for(size_t i=0; i<m_npop; i++) {
        m_config.push_back(config_t());
        m_logbf.push_back(0.0);
    }
   
    // initialize mvb log likelihood
    double denom = log_mvb_denom(params);
    m_loglik_mvb = m_params(0)*((double)m_nsnp)-denom*(double(m_nsnp));
}

// update the state
void Chain::next_state() {
   
    // choose an index to flip at random
    size_t idx = rand_int(0, m_nsnp-1);
    if(m_nsnp == 0) idx = 0;

    // get next possible states
    config_t pop1_next0 = m_config[0]; pop1_next0.erase(idx);
    config_t pop1_next1 = m_config[0]; pop1_next1.insert(idx);
    config_t pop2_next0 = m_config[1]; pop2_next0.erase(idx);
    config_t pop2_next1 = m_config[1]; pop2_next1.insert(idx);

    vector<config_t> next_conf00, next_conf01, next_conf10, next_conf11;
    next_conf00.push_back(pop1_next0); next_conf00.push_back(pop2_next0);
    next_conf01.push_back(pop1_next1); next_conf01.push_back(pop2_next0);
    next_conf10.push_back(pop1_next0); next_conf10.push_back(pop2_next1);
    next_conf11.push_back(pop1_next1); next_conf11.push_back(pop2_next1);

    // evaluate the likelihood of the next states
    double llk_mvb00 = eval_loglik_mvb(next_conf00);
    double llk_mvb01 = eval_loglik_mvb(next_conf01);
    double llk_mvb10 = eval_loglik_mvb(next_conf10);
    double llk_mvb11 = eval_loglik_mvb(next_conf11);

    double logbf_pop1_next0;
    if(pop1_next0 == m_config[0]) logbf_pop1_next0 = m_logbf[0];
    else                          logbf_pop1_next0 = log_bf(0, pop1_next0);
    
    double logbf_pop1_next1;
    if(pop1_next1 == m_config[0]) logbf_pop1_next1 = m_logbf[0];
    else                          logbf_pop1_next1 = log_bf(0, pop1_next1);
    
    double logbf_pop2_next0;
    if(pop2_next0 == m_config[1]) logbf_pop2_next0 = m_logbf[1];
    else                          logbf_pop2_next0 = log_bf(1, pop2_next0);
    
    double logbf_pop2_next1;
    if(pop2_next1 == m_config[1]) logbf_pop2_next1 = m_logbf[1];
    else                          logbf_pop2_next1 = log_bf(1, pop2_next1);

    // get log likelihood for next four states
    double llk_zsc00 = logbf_pop1_next0 + logbf_pop2_next0 + m_lambda;
    double llk_zsc01 = logbf_pop1_next1 + logbf_pop2_next0;
    double llk_zsc10 = logbf_pop1_next0 + logbf_pop2_next1;
    double llk_zsc11 = logbf_pop1_next1 + logbf_pop2_next1;

    // incorporate temperature
    double llk00 = (llk_mvb00 + llk_zsc00)/m_temp;
    double llk01 = (llk_mvb01 + llk_zsc01)/m_temp;
    double llk10 = (llk_mvb10 + llk_zsc10)/m_temp;
    double llk11 = (llk_mvb11 + llk_zsc11)/m_temp;

    // convert to probability
    double prob_denom = logsum(logsum(logsum(llk00, llk01), llk10), llk11);
    double prob00 = exp(llk00 - prob_denom);
    double prob01 = exp(logsum(llk00, llk01) - prob_denom);
    double prob10 = exp(logsum(llk10, logsum(llk00, llk01)) - prob_denom);
    
    // advance to next state
    double rnd = rand_double();
    if(rnd < prob00) {
        m_config = next_conf00;
        m_logbf[0] = logbf_pop1_next0;
        m_logbf[1] = logbf_pop2_next0;
    }
    else if(rnd >= prob00 && rnd < prob01) {
        m_config = next_conf01;
        m_logbf[0] = logbf_pop1_next1;
        m_logbf[1] = logbf_pop2_next0;
    }
    else if(rnd >= prob01 && rnd < prob10) {
        m_config = next_conf10;
        m_logbf[0] = logbf_pop1_next0;
        m_logbf[1] = logbf_pop2_next1;
    }
    else {
        m_config = next_conf11;
        m_logbf[0] = logbf_pop1_next1;
        m_logbf[1] = logbf_pop2_next1;
    }

    // increment counter
    m_nstep++; 
}

// compute the log of bayes factor
double Chain::log_bf(size_t updt, const config_t& conf) {

    // initialization
    size_t num_cau = conf.size();
    if(num_cau == 0) return 0.0;

    // extract the z-score and ld
    VectorXf zsc_conf(num_cau);
    MatrixXf ld_conf = MatrixXf::Zero(num_cau, num_cau);
    size_t ii = 0;
    for (auto it1=conf.begin(); it1 != conf.end(); it1++) {
        size_t i = *it1;
        zsc_conf(ii) = m_zsc[updt](i);
        size_t jj = 0;
        for (auto it2=conf.begin(); it2 != conf.end(); it2++) {
            size_t j = *it2;
            ld_conf(jj, ii) = m_ld[updt](j, i);
            jj++;
        }
        ii++;
    }

    // compute svd
    BDCSVD<MatrixXf> svd(ld_conf, ComputeFullV);
    MatrixXf D = svd.singularValues().asDiagonal();
    MatrixXf V = svd.matrixV();

    // find rank of the matrix
    size_t num_indep = 0.0;
    for(size_t i=0; i<num_cau; i++)
        if(D(i,i) > 0.0) num_indep++; else break;

    // compute log Bayes factor
    double sigma_sq = m_sigma_sq[updt]/((double) num_indep);
    double bf = 0.0;
    for(size_t i=0; i<num_indep; i++) {
        bf -= log(1.0 + sigma_sq*D(i,i));
        double prod = zsc_conf.transpose()*V.col(i);
        bf += sigma_sq/(1.0+sigma_sq*D(i,i))*(prod*prod);
    }
    bf *= 0.5;

    return bf;
}

// evaluate the loglikelihood of the configuration
double Chain::eval_loglik_mvb(vector<config_t> config) {
    
    config_t conf1 = config[0];
    config_t conf2 = config[1];

    double denom = log_mvb_denom(m_params);
    double llk = 0.0 - denom*((double)m_nsnp);

    double n11 = 0.0;
    for(auto it=conf1.begin(); it != conf1.end(); it++) {
        n11 += (conf2.find(*it) != conf2.end());
    }
    double n01 = conf1.size() - n11;
    double n10 = conf2.size() - n11;
    double n00 = m_nsnp - n01 - n10 - n11;
    llk = llk + n00*(m_params(0))+
                n01*(m_params(1)+m_params(0))+
                n10*(m_params(2)+m_params(0))+
                n11*(m_params(3)+m_params(2)+m_params(1)+m_params(0));
    
    return llk;
}

// get state
vector<config_t> Chain::get_state() {
    return m_config;
}

// set state
void Chain::set_state(const vector<config_t>& state) {
    for(size_t i=0; i<state.size(); i++) {
        m_config[i] = state[i];
    }
}

// return the mvb log likelihood
double Chain::get_loglik_mvb() {
    return m_loglik_mvb;
}

// return the temperature
double Chain::get_temp() {
    return m_temp;
}

// reset the chain
void Chain::reset(const VectorXf& params) {
    
    // clear the counter
    m_nstep = 0;
    m_config.clear();
    m_logbf.clear();

    // initialize the configs, mean, and loglik for each zscore
    for(size_t i=0; i<m_npop; i++) {
        m_config.push_back(config_t());
        m_logbf.push_back(0.0);
    }
    
    // initialize mvb log likelihood
    m_params = params;
    double denom = log_mvb_denom(params);
    m_loglik_mvb = m_nsnp*m_params(0) - denom*((double)m_nsnp);
}

// create a sampler and burn
Sampler::Sampler(const vector<Chain>& chains, size_t nburn, size_t nsample):
    m_chains(chains), m_nburn(nburn), m_nsample(nsample), m_ndone(0) {

    size_t nchain = m_chains.size();
    for(size_t i=0; i<nburn; i++) {
        for(size_t j=0; j < nchain; j++) {
            m_chains[j].next_state();
        }
        swap_chain();
    }
}

// sample from the sampler
vector<config_t> Sampler::sample() {
 
    // advance chain
    size_t nchain = m_chains.size();
    for(size_t i=0; i < nchain; i++)
        m_chains[i].next_state();
    swap_chain();
    m_ndone++;
    
    return m_chains[0].get_state();
}

// reset the sampler with new params
void Sampler::reset(const VectorXf& params) {
    
    // reset chain
    m_ndone = 0;
    for(size_t i=0; i<m_chains.size(); i++) {
        m_chains[i].reset(params);
    }
    
    // burn in 
    size_t nchain = m_chains.size();
    for(size_t i=0; i<m_nburn; i++) {
        for(size_t j=0; j < nchain; j++) {
            m_chains[j].next_state();
        }
        swap_chain();
    }
}

// check if need to generate more sample
bool Sampler::done() {
    return m_ndone == m_nsample;
}

// swap the state of two chain
void Sampler::swap_chain() {
    
    size_t nchain = m_chains.size();
    if(nchain < 2) return;
    size_t i = m_ndone % (nchain-1);
    size_t j = i + 1;
    vector<config_t> state_i = m_chains[i].get_state();
    vector<config_t> state_j = m_chains[j].get_state();

    double pii = m_chains[i].get_loglik_mvb();
    double pjj = m_chains[j].get_loglik_mvb();
    double pij = m_chains[i].eval_loglik_mvb(state_j);
    double pji = m_chains[j].eval_loglik_mvb(state_i);
            
    double prob = (pij+pji)-(pii+pjj);
    prob = prob > 0 ? 1 : exp(prob);
    if(rand_double() < prob) {
        vector<config_t> tmp = state_i;
        m_chains[i].set_state(state_j);
        m_chains[j].set_state(tmp);
    }
}

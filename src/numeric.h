#include <math.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm>

#include "utils.h"

#ifndef NUMERIC
#define NUMERIC

typedef set<size_t> config_t;

// define epsilon
const double eps = 1e-8;

// define sub set
const size_t sub0[] = {0};
const vector<size_t> sub0_vec(sub0, sub0+sizeof(sub0)/sizeof(size_t)); 
const size_t sub1[] = {0, 1};
const vector<size_t> sub1_vec(sub1, sub1+sizeof(sub1)/sizeof(size_t));
const size_t sub2[] = {0, 2};
const vector<size_t> sub2_vec(sub2, sub2+sizeof(sub2)/sizeof(size_t));
const size_t sub3[] = {0, 1, 2, 3};
const vector<size_t> sub3_vec(sub3, sub3+sizeof(sub3)/sizeof(size_t));
const vector<size_t> subset_arr[]={sub0_vec,sub1_vec,sub2_vec,sub3_vec};
const vector<vector<size_t> > subset(subset_arr,
    subset_arr+sizeof(subset_arr)/sizeof(vector<size_t>));

// define sup set
const size_t sup0[] = {0, 1, 2, 3};
const vector<size_t> sup0_vec(sup0, sup0+sizeof(sup0)/sizeof(size_t)); 
const size_t sup1[] = {1, 3};
const vector<size_t> sup1_vec(sup1, sup1+sizeof(sup1)/sizeof(size_t));
const size_t sup2[] = {2, 3};
const vector<size_t> sup2_vec(sup2, sup2+sizeof(sup2)/sizeof(size_t));
const size_t sup3[] = {3};
const vector<size_t> sup3_vec(sup3, sup3+sizeof(sup3)/sizeof(size_t));
const vector<size_t> supset_arr[]={sup0_vec,sup1_vec,sup2_vec,sup3_vec};
const vector<vector<size_t> > supset(supset_arr,
    supset_arr+sizeof(supset_arr)/sizeof(vector<size_t>));

// compute log density of a multivariate normal
double mvn_logpdf(const VectorXf& mean, const MatrixXf& sigma,
    const VectorXf& sample);

// compute the correlation between two vectors
double corr(const doublevec_t& v1, const doublevec_t& v2);

// compute the ld matrix from the genotye data
MatrixXf compute_ld(const doublemat_t& gen_mat);

// compute log som
double logsum(double a, double b);

// mvb s function
double mvb_sfunc(const VectorXf& params, size_t c);

// log mvb denom function
double log_mvb_denom(const VectorXf& params);

// return true of a is a subset of b
bool is_asubb(size_t a, size_t b);

// estimate causal effect
VectorXf estimate_caueff(const MatrixXf& ld_mat, const VectorXf& zsc,
    double lambda, VectorXf& caueff_var, double& hsq);

// estimate sigmasq 
double estimate_sigmasq(const MatrixXf& ld_mat, const VectorXf& zsc,
    double lambda);

#endif

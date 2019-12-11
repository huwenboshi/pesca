#include <vector>
#include "utils.h"
#include <iostream>
#include <fstream>

#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif
#include <Eigen/Dense>
#include "mcmc.h"

using namespace std;

#ifndef FIT
#define FIT

// update parameters
void update_params(VectorXf& params, const VectorXf& conf_cnt);

// fit the model
void fit(const vector<vector<VectorXf> >& all_zsc,
         const vector<vector<double> >& all_sigma_sq,
         const vector<vector<MatrixXf> >& all_ld,
         const vector<size_t>& all_nsnp, VectorXf& params,
         double templ, double temph, size_t nchain,
         size_t nburn, size_t nsample, double lambda,
         size_t max_iter, ofstream& outfile, ofstream& logfile);

#endif

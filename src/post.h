#include <vector>
#include "utils.h"
#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif
#include <Eigen/Dense>
#include "mcmc.h"

using namespace std;

#ifndef POST
#define POST

// fit the model
void post(const vector<vector<VectorXf> >& all_zsc,
         const vector<vector<double> >& all_sigma_sq,
         const vector<vector<MatrixXf> >& all_ld,
         const vector<size_t>& all_nsnp, VectorXf& params,
         double templ, double temph, size_t nchain,
         size_t nburn, size_t nsample, double lambda, size_t min_iter, 
         const vector<vector<Zscore_t> >& all_zsc_inf1,
         const vector<vector<Zscore_t> >& all_zsc_inf2,
         ofstream& outfile, ofstream& logfile);

#endif

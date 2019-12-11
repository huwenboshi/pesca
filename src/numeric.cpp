#include <time.h>
#include "numeric.h"

using namespace std;

// compute log density of a multivariate normal
double mvn_logpdf(const VectorXf& mean, const MatrixXf& sigma,
    const VectorXf& sample) {
    
    double logpdf = 0.0;
    double dim = sample.rows();
    MatrixXf sigma_inv = sigma.ldlt().solve(MatrixXf::Identity(dim, dim));
    VectorXf diff = sample - mean;
    logpdf = -0.5 * diff.transpose() * sigma_inv  * diff;
    
    double dtm = sigma.determinant();
    if(dtm < 10e-8) return -10e8;
    
    logpdf = logpdf - 0.5 * log(dtm);
    logpdf = logpdf - 0.5 * dim * log(2.0*M_PI);
    
    return logpdf;
}

// compute the correlation between two vectors
double corr(const doublevec_t& v1, const doublevec_t& v2) {

    double mean1=0.0, mean1_sq=0.0;
    double mean2=0.0, mean2_sq=0.0;
    double mean12 = 0.0;
    size_t num_indv = v1.size();

    for(size_t i=0; i < num_indv; i++) {
        mean1 += v1[i];
        mean1_sq += v1[i] * v1[i];
        mean2 += v2[i];
        mean2_sq += v2[i] * v2[i];
        mean12 += v1[i] * v2[i];
    }

    mean1 /= (double) num_indv;
    mean1_sq /= (double) num_indv;
    mean2 /= (double) num_indv;
    mean2_sq /= (double) num_indv;
    mean12 /= (double) num_indv;

    double cov = mean12 - mean1*mean2;
    double var1 = mean1_sq - mean1*mean1;
    double var2 = mean2_sq - mean2*mean2;

    if(var1*var2 == 0.0) {
        return 0.0;
    }

    double r = cov / sqrt(var1 * var2);
    return r;
}

// compute the ld matrix from the genotye data
MatrixXf compute_ld(const doublemat_t& gen_mat) {

    size_t num_snps = gen_mat.size();
    MatrixXf ld_mat(num_snps, num_snps);

    for(size_t i=0; i < num_snps; i++) {
        ld_mat(i,i) = 1.0;
        for(size_t j=i+1; j < num_snps; j++) {
            double r = corr(gen_mat[i], gen_mat[j]);
            ld_mat(i,j) = r;
            ld_mat(j,i) = r;
        }
    }

    return ld_mat;
}

// compute log som
double logsum(double a, double b) {

    if(a > b) {
        double c = a;
        a = b;
        b = c;
    }

    return b + log(exp(a-b) + 1.0);
}

// mvb s function
double mvb_sfunc(const VectorXf& params, size_t c) {

    double sum = 0.0;

    for(size_t i=0; i < subset[c].size(); i++) {
        size_t idx = subset[c][i];
        sum += params(idx);
    }

    return sum;
}

// log mvb denom function
double log_mvb_denom(const VectorXf& params) {

    double sum = 0.0;

    for(size_t i=0; i < params.rows(); i++) {
        if(i == 0) {
            sum = mvb_sfunc(params, i);
        }
        else {
            sum = logsum(sum, mvb_sfunc(params, i));
        }
    }

    return sum;
}

// check if a is a subset of b
bool is_asubb(size_t a, size_t b) {
    
    if(a == 0) {
        return true;
    }

    if(a == 1) {
        if(b == 3 || b == 1) {
            return true;
        }
    }

    if(a == 2) {
        if(b == 3 || b == 2) {
            return true;
        }
    }

    if(a == 3) {
        if(b == 3) {
            return true;
        }
    }

    return false;
}

// estimate causal effect scaled by n
VectorXf estimate_caueff(const MatrixXf& ld_mat, const VectorXf& zsc, 
    double lambda, VectorXf& caueff_var, double& sigmasq) {
    
    size_t nsnp = zsc.rows();
   
    BDCSVD<MatrixXf> svd(ld_mat, ComputeFullU | ComputeFullV);
    MatrixXf D = svd.singularValues().asDiagonal();
    MatrixXf U = svd.matrixU();
    MatrixXf V = svd.matrixV();

    MatrixXf mat_inv(ld_mat.rows(), ld_mat.cols());
    mat_inv.setZero();
    for(int i=0; i < nsnp; i++) {
        if(D(i,i) < lambda) continue;
        mat_inv += 1.0/D(i,i) * U.col(i) * V.col(i).transpose();
    }

    caueff_var = mat_inv.diagonal();
    VectorXf caueff = mat_inv*zsc;
    sigmasq = 0.0;
    for(int i=0; i<caueff.rows(); i++) {
        sigmasq += caueff(i)*caueff(i);
    }
    return caueff;
}

// estimate sigmasq
double estimate_sigmasq(const MatrixXf& ld_mat, const VectorXf& zsc, 
    double lambda) {
    
    size_t nsnp = zsc.rows();
   
    BDCSVD<MatrixXf> svd(ld_mat, ComputeFullU | ComputeFullV);
    MatrixXf D = svd.singularValues().asDiagonal();
    MatrixXf U = svd.matrixU();
    MatrixXf V = svd.matrixV();

    MatrixXf mat_inv(ld_mat.rows(), ld_mat.cols());
    mat_inv.setZero();
    for(int i=0; i < nsnp; i++) {
        if(D(i,i) < lambda) continue;
        mat_inv += 1.0/D(i,i) * U.col(i) * V.col(i).transpose();
    }

    VectorXf caueff = mat_inv*zsc;
    double sigmasq = 0.0;
    for(int i=0; i<caueff.rows(); i++) {
        sigmasq += caueff(i)*caueff(i);
    }
    
    return sigmasq;
}

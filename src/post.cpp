#include "post.h"
#include "numeric.h"

using namespace std;

#define DEBUG_FIT true

void post(const vector<vector<VectorXf> >& all_zsc,
         const vector<vector<double> >& all_sigma_sq,
         const vector<vector<MatrixXf> >& all_ld,
         const vector<size_t>& all_nsnp, VectorXf& params,
         double templ, double temph, size_t nchain,
         size_t nburn, size_t nsample, double lambda, size_t min_iter,
         const vector<vector<Zscore_t> >& all_zsc_inf1,
         const vector<vector<Zscore_t> >& all_zsc_inf2,
         ofstream& outfile, ofstream& logfile) {

    // count total number of SNPs in the loci
    double nsnp_tot = 0.0;
    for(size_t l=0; l<all_nsnp.size(); l++) {
        nsnp_tot += all_nsnp[l];
    }
    size_t nloci = all_zsc.size();

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

    // iterate through loci
    outfile << "SNP\tBP\tA0_GWAS1\tA1_GWAS1\tZ_GWAS1\tA0_GWAS2\t"
         << "A1_GWAS2\tZ_GWAS2\tP_GWAS1_ONLY\tP_GWAS2_ONLY\tP_BOTH" << endl;
    for(size_t l=0; l<nloci; l++) {
        
        // record the number of configurations
        double* conf_cnt = new double[all_nsnp[l]*3];
        for(size_t idx=0; idx<all_nsnp[l]*3; idx++)
            conf_cnt[idx] = 0.0;
        for(size_t iter=0; iter < min_iter; iter++) {
           
            // get samples
            while(!all_sampler[l].done()) {

                // get a sample
                vector<config_t> sample = all_sampler[l].sample();
                config_t conf1 = sample[0];
                config_t conf2 = sample[1];

                // count the occurance
                for(auto it=conf1.begin(); it != conf1.end(); it++) {
                    conf_cnt[(*it)*3+0] += 1.0;
                }
                for(auto it=conf2.begin(); it != conf2.end(); it++) {
                    conf_cnt[(*it)*3+1] += 1.0;
                }
                for(auto it=conf1.begin(); it != conf1.end(); it++) {
                    if(conf2.find(*it) != conf2.end()) {
                        conf_cnt[(*it)*3+2] += 1.0;
                    }
                }
            }
            // reset chain for next round
            all_sampler[l].reset(params);
        }

        // get posterior
        for(size_t i = 0; i < all_nsnp[l]; i++) {
            double p1 = conf_cnt[i*3+0] / nsample / min_iter;
            double p2 = conf_cnt[i*3+1] / nsample / min_iter;
            double p12 = conf_cnt[i*3+2] / nsample / min_iter;
            outfile << all_zsc_inf1[l][i].m_rsid << "\t"
                 << all_zsc_inf1[l][i].m_pos  << "\t"
                 << all_zsc_inf1[l][i].m_ref  << "\t"
                 << all_zsc_inf1[l][i].m_alt  << "\t"
                 << all_zsc_inf1[l][i].m_zsc  << "\t"
                 << all_zsc_inf2[l][i].m_ref  << "\t"
                 << all_zsc_inf2[l][i].m_alt  << "\t"
                 << all_zsc_inf2[l][i].m_zsc  << "\t"
                 << p1-p12 << "\t" << p2-p12 << "\t" << p12 << endl;
        }

        delete[] conf_cnt;
    }
}

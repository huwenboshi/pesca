#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif
#include<Eigen/Dense>
#include <boost/program_options.hpp>

using namespace Eigen;
using namespace std;

namespace po = boost::program_options;

#ifndef UTILS
#define UTILS

typedef vector<string> strvec_t;
typedef vector<double> doublevec_t;
typedef vector<vector<double> > doublemat_t;

typedef struct Zscore_t {
    string m_rsid;
    size_t m_pos;
    char m_ref;
    char m_alt;
    double m_zsc;
    double m_n;
} Zscore_t;

// split string
strvec_t split_string(const string& str);

// get command line input
po::variables_map get_command_line(int ac, char* av[]);

// load zscore file
vector<Zscore_t> load_zscore(const char* zsc_file_nm);

// load zscore list
vector<string> load_list(const char* list_nm);

// load genotype file
doublemat_t load_genotype(const char* gen_file_nm);

MatrixXf load_ldmat(const char* ldmat_file_nm);

// doublevec to vectorxd
VectorXf get_zsc_vec(const vector<Zscore_t>& v);



#endif

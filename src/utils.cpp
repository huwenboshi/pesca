#include "utils.h"

using namespace std;

// split string
vector<string> split_string(const string& str) {
   
    vector<string> split;
    stringstream ssin(str);
    
    while(ssin.good()) {
        string tmp;
        ssin >> tmp;
        split.push_back(tmp);
    }

    return split;
}

// get command line input
po::variables_map get_command_line(int ac, char* av[]) {
    
    po::options_description desc("Allowed options");
    
    desc.add_options()
        ("help", "produce help message")
        ("zscore1", po::value<string>(), "summary stats file 1")
        ("zscore2", po::value<string>(), "summary stats file 2")
        ("ld1", po::value<string>(), "ld 1")
        ("ld2", po::value<string>(), "ld 2")
        ("mode", po::value<string>(), "mode")
        ("sigmasq1", po::value<double>(), "sigmasq1")
        ("sigmasq2", po::value<double>(), "sigmasq2")
        ("totnsnp", po::value<double>(), "totnsnp")
        ("templ", po::value<double>()->default_value(1.0),
            "lower temperature for mcmc")
        ("temph", po::value<double>()->default_value(1.0),
            "higher temperature for mcmc")
        ("nchain", po::value<size_t>()->default_value(1),
            "number of markov chains")
        ("nburn", po::value<size_t>()->default_value(30000),
            "number of burn in")
        ("nsample", po::value<size_t>()->default_value(50000),
            "number of samples")
        ("f00", po::value<double>()->default_value(0.0), "f00")
        ("f01", po::value<double>()->default_value(-3.9), "f01")
        ("f10", po::value<double>()->default_value(-3.9), "f10")
        ("f11", po::value<double>()->default_value(3.9), "f11")
        ("lambda", po::value<double>()->default_value(0.0001), "lambda")
        ("max_iter", po::value<size_t>()->default_value(100), "max_iter")
        ("out", po::value<string>(), "out")
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line(ac, av, desc), vm);
    try {
        po::notify(vm);
    } catch (exception& e) {
        cerr << "test" << endl;
        cerr << desc << endl;
        exit(1);
    }
    return vm;
}

// load zscores
vector<Zscore_t> load_zscore(const char* zsc_file_nm) {
    
    vector<Zscore_t> zscore_info;

    ifstream zsc_file(zsc_file_nm);
    if(zsc_file.fail()){
        cerr << "Error: Failed to open " << string(zsc_file_nm) << "!" << endl;
        exit(1);
    }
    string line;
    getline(zsc_file, line);
    while(getline(zsc_file, line)) {
        strvec_t split = split_string(line);
        Zscore_t tmp;
        tmp.m_rsid = split[0];
        tmp.m_pos = atol(split[1].c_str());
        tmp.m_ref = split[2][0];
        tmp.m_alt = split[3][0];
        tmp.m_zsc = atof(split[4].c_str());
        zscore_info.push_back(tmp);
    }
    zsc_file.close();

    return zscore_info;
}

// load zscore list
vector<string> load_list(const char* list_nm) {
    
    vector<string> list;

    ifstream list_file(list_nm);
    if(list_file.fail()){
        cerr << "Error: Failed to open " << string(list_nm) << "!" << endl;
        exit(1);
    }
    string line;
    while(getline(list_file, line)) {
        list.push_back(line);
    }
    list_file.close();

    return list;

}

// load genotype file
doublemat_t load_genotype(const char* gen_file_nm) {

    vector<vector<double> > gen_mat;

    ifstream gen_file(gen_file_nm);
    if(gen_file.fail()){
        cerr << "Error: Failed to open " << string(gen_file_nm) << "!" << endl;
        exit(1);
    }
    string line;
    while(getline(gen_file, line)) {
        strvec_t split = split_string(line);
        doublevec_t tmp;
        for(size_t i=0; i < split.size(); i++) {
            tmp.push_back(atof(split[i].c_str()));
        }
        gen_mat.push_back(tmp);
    }
    gen_file.close();

    return gen_mat;
}

// load genotype file
MatrixXf load_ldmat(const char* ldmat_file_nm) {

    vector<vector<double> > ldmat;

    ifstream ldmat_file(ldmat_file_nm);
    if(ldmat_file.fail()){
        cerr << "Error: Failed to open " << string(ldmat_file_nm)<<"!" << endl;
        exit(1);
    }
    string line;
    while(getline(ldmat_file, line)) {
        strvec_t split = split_string(line);
        doublevec_t tmp;
        for(size_t i=0; i < split.size(); i++) {
            tmp.push_back(atof(split[i].c_str()));
        }
        ldmat.push_back(tmp);
    }
    ldmat_file.close();

    MatrixXf ldmat_mat(ldmat.size(), ldmat.size());
    for(size_t i=0; i<ldmat_mat.rows(); i++) {
        for(size_t j=0; j<ldmat_mat.cols(); j++) {
            ldmat_mat(i,j) = ldmat[i][j];
        }
    }

    return ldmat_mat;
}

// doublevec to vectorxd
VectorXf get_zsc_vec(const vector<Zscore_t>& v) {

    size_t num_snps = v.size();
    VectorXf vec(num_snps);

    for(size_t i=0; i < num_snps; i++) {
        vec(i) = v[i].m_zsc;
    }

    return vec;
}

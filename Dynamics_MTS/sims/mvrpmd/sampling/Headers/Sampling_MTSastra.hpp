#ifndef Sampling_MTSastra_hpp
#define Sampling_MTSastra_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include "MVRPMD_MTS_Hamiltonian.hpp"

#include "nr3.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"

using namespace boost::numeric::ublas;

class Sampling_MTSastra{
    
public:
    
    Sampling_MTSastra(int my_id, int num_procs, int root_proc, int num_beads, int elec_beads, double mass,
                 int num_states, double beta_num_beads, double beta_elec_beads, double nuc_ss,
                 double elec_ss, std::string root);
    
    /* Run Sampling simulation. A batch of samples will be produces after
     running the function.*/
    void runSimulation();
    
    void set_decor_len(unsigned long long set_In);
    
    void set_num_samples(int set_In);
    
private:
    
    /* Private functions. */
    void sample_nuc();
    
    void sample_elec();
    
    inline double nuc_dist(const double rn);
    
    inline double elec_dist(const double rn);
    
    inline int rand_nuc_bead(const Ullong rn);
    
    inline int rand_elec_bead(const Ullong rn);
    
    void readPSV();
    
    void stash_Q(int sample);
    
    void stash_x(int sample);

    void stash_p(int sample);

    void print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt);

    void print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt);

    void saveQ();

    void saveP();

    void savex();

    void savep();

    template <typename T>
    void histogram(std::string fileName, int bins, T X);

    /* Private Data */
    int my_id, num_procs, root_proc;

    double mass;
    int num_beads; //number of nuclear beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    double beta_num_beads; //beta/num_beads
    double energy; //energy after a proposed move
    double nuc_ss; //nuclear step size
    double elec_ss; //electronc step size

    unsigned long long decor_len; //decorrelation length
    int num_samples; //number of samples

    unsigned long long sys_steps_accpt; //system steps accepted
    unsigned long long sys_steps; //system steps accepted

    unsigned long long elec_steps_accpt; //system steps accepted
    unsigned long long elec_steps; //system steps accepted

    vector<double> Q; //vector of bead positions
    vector<double> Q_prop; //vector of bead positions
    vector<double> Q_samp; //vector of all samples

    /* matrix of electronic x variables; x[i][j] ith bead and jth state*/
    matrix<double> x;
    matrix<double> x_prop; //proposed move of x variables
    matrix<double> x_samp; //matrix of all samples

    /* matrix of electronic p variables; p[i][j] ith bead and jth state*/
    matrix<double> p;
    matrix<double> p_prop; //proposed move of p variables
    matrix<double> p_samp; //matrix of all samples

    std::string root;

    SpringEnergy V_spring;
    StateIndepPot V0;
    GTerm G;
    M_Matrix M;
    M_Matrix_MTS M_MTS;
    C_Matrix C;
    Theta_MTS thetaMTS;
    MVRPMD_MTS_Hamiltonian H_MTS;
    
    Ran myRand; //NR3 random number generator
    
};

#endif

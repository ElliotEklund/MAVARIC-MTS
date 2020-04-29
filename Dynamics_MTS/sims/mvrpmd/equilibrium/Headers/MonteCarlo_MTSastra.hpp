#ifndef MonteCarlo_MTSastra_hpp
#define MonteCarlo_MTSastra_hpp

#include <stdio.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <string>

#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "MVRPMD_MTS_Estimator.hpp"
#include "MonteCarloHelper.h"

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class MonteCarlo_MTSastra {
    
public:
    MonteCarlo_MTSastra(int my_id, int root_proc, int num_procs, int num_beads, int elec_beads, double mass,
                   int num_states,double beta_num_beads, double beta_elec_beads, double nuc_ss,
                   double elec_ss,std::string root);

    /* Generate initial conditions for Q. Q is set after calling*/
    void gen_initQ();

    /* Generate initial conditions for x. x is set after calling*/
    void gen_initx();

    /* Generate initial conditions for p. p is set after calling*/
    void gen_initp();

    /* Run monte carlo simulation.*/
    void runSimulation();

    /* Mutators: Function is self-evident */
    void set_num_steps(unsigned long long num_steps_In);
    
    void set_esti_rate(unsigned long long esti_rate_In);
    
    void set_write_PSV(bool set_In);
    
    void set_read_PSV(bool set_In);
    
    void set_read_Data(bool set_In);

    void set_write_Data(bool set_In);
    
    void print_avg_energy(double estimator_total, double sgn_total);

private:
    
    /* Private functions. */
    void sample_nuc();
    void sample_elec();

    inline double nuc_dist(const double rn);
   
    inline double elec_dist(const double rn);
  
    inline int rand_nuc_bead(const Ullong rn);
   
    inline int rand_elec_bead(const Ullong rn);

    SpringEnergy V_spring;
    StateIndepPot V0;
    GTerm G;
    M_Matrix M;
    M_Matrix_MTS M_MTS;
    C_Matrix C;
    Theta_MTS thetaMTS;
    dTheta_MTS_dBeta dthetaMTS_dBeta;
    MVRPMD_MTS_Hamiltonian H_MTS;
    MVRPMD_MTS_Estimator Esti_MTS;

    /* Private data. */
    int my_id; //unique processor id
    int root_proc; //root processor
    int num_procs; //number of processors
    
    bool writePSV; //write PSV to file after simulation if true
    bool readPSV; //read in PSC.txt before simulation if true
    bool readData; //read in mcData.txt before simulation if true
    bool writeData; //write to mcData.txt after simulation if true
    
    int num_beads; //number of nuclear beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    double beta_num_beads; //beta/num_beads
    double energy; //energy after a proposed move
    double estimator; //energy estimator after a proposed move
    double sgn_theta; //sign of theta
    double nuc_ss; //nuclear step size
    double elec_ss; //electronic step size

    unsigned long long num_steps; //number of monte carlo steps
    unsigned long long esti_rate; //rate at which estimator is calculated

    unsigned long long sys_steps_accpt; //system steps accepted
    unsigned long long sys_steps; //system steps accepted

    unsigned long long elec_steps_accpt; //system steps accepted
    unsigned long long elec_steps; //system steps accepted

    vector<double> Q; //vector of bead positions
    vector<double> Q_prop; //vector of bead positions


    /* matrix of electronic x variables; x[i][j] ith bead and jth state*/
    matrix<double> x;
    matrix<double> x_prop; //proposed move of x variables

    /* matrix of electronic p variables; p[i][j] ith bead and jth state*/
    matrix<double> p;
    matrix<double> p_prop; //proposed move of p variables

    vector<double> estimator_t; //estimator evaluated at "time t"

    Ran myRand; //NR3 random number generator


    MonteCarloHelper myHelper;
};

#endif


#ifndef equilib_mvrpmd_hpp
#define equilib_mvrpmd_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <stdlib.h>
#include <time.h>
#include <string>

#include "nr3.h"
#include "ran.h"
#include "mpi.h"

#include "mpi_wrapper.hpp"
#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "mvrpmd_mixed.hpp"
#include "mvrpmd_mixed_ham.hpp"
#include "MVRPMD_MTS_Estimator.hpp"
#include "mvrpmd_mixed_esti.hpp"
#include "MonteCarloHelper.h"
#include "system_step.hpp"
#include "elec_step.hpp"

class equilib_mvrpmd{
    
public:
    equilib_mvrpmd(int my_id, int root_proc, int num_procs, std::string root_path);
    
    /* Runs Monte Carlo simulation.
     nuc_ss: nuclear step size used by nuclear part of system; used for proposing
              monte carlo steps
     x_ss: x mapping variable step size used by x part of system; used for proposing
              monte carlo steps
     p_ss: p mapping variable step size used by p part of system; used for proposing
              monte carlo steps
     num_steps: number of monte carlo steps simulation will take
     stride: length between successive collection of energy estimator
     */
    void run(double nuc_ss, double x_ss, double p_ss,unsigned long long num_steps,
             unsigned long long stride);

/* Mutators */
    void initialize_system(int nuc_beads_IN,int elec_beadsIN, int num_statesIN,
                           double massIN,double betaIN,double alpha);
    
    void initialize_files(bool writePSV_IN, bool readPSV_IN,
                          bool writeData_IN, bool readData_IN);

    void set_write_PSV(bool set_In);

    void set_read_PSV(bool set_In);

    void set_read_Data(bool set_In);

    void set_write_Data(bool set_In);

private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    bool writePSV, readPSV; // determine how Phase Spave Variables are processed
    bool writeData, readData; // determine how MC data is processed
    
    /* Number of nuclear beads, electronic beads, and electronic states*/
    int nuc_beads, elec_beads, num_states;
    double mass, beta; //system mass; 1/kb T
    double alpha; //mapping variable prefactor
    
    bool sys_set, files_set; //true only once initialize_system(initialize_files) called

/* Objects */
    Ran myRand; //NR3 random number generator
    
    MonteCarloHelper helper;//helper object

/* Functions */
    /* Initialize Q with random vaules.
     Q: vector of bead positions
     num_beads: number of beads in Q
     step_size: sets range of individual elements of Q; Q[i] ~ [-step_size,step_size]*/
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);
    
    /* Initialize v with random vaules.
     v: matrix of mapping variables (either x or p)
     num_beads: number of beads in v
     num_states: number of electronic states in v
     step_size: sets range of individual elements of v; v[i,j] ~ [-step_size,step_size]*/
    void gen_initElec(matrix<double> &v, int num_beads, int num_states,
                      double step_size);
    
    /* Returns a uniform random number between [-step_size,step_size]
     rn: random double
     step_size: specifices range of uniform random distribution*/
    inline double step_dist(const double rn, double step_size);
};

#endif

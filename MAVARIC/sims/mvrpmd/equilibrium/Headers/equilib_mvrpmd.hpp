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
#include "MVRPMD_MTS_Estimator.hpp"
#include "MonteCarloHelper.h"
#include "system_step.hpp"
#include "elec_step.hpp"

class equilib_mvrpmd{
    
public:
    equilib_mvrpmd(int my_id, int root_proc, int num_procs, std::string root_path);
    
    void run(double nuc_ss, double x_ss, double p_ss,unsigned long long num_steps,
             unsigned long long stride);

/* Mutators */
    void initialize_system(int nuc_beads_IN,int elec_beadsIN, int num_statesIN,
                           double massIN,double betaIN);
    
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
        
    int nuc_beads, elec_beads, num_states;
    double mass, beta;

/* Objects */
    Ran myRand; //NR3 random number generator
    
    MonteCarloHelper helper;

/* Functions */
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);
    
    void gen_initElec(matrix<double> &v, int num_beads, int num_states,
                      double step_size);
    
    inline double step_dist(const double rn, double step_size);
    

};

#endif

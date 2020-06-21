#ifndef dynamics_mvrpmd_hpp
#define dynamics_mvrpmd_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "auto_correlation.hpp"
#include "init_PAC.hpp"
#include "energy_conserv.hpp"

#include <string>
#include "mpi.h"

using namespace boost::numeric::ublas;

class dynamics_mvrpmd{
    
public:
    
    dynamics_mvrpmd(int my_id, int num_procs, int root_proc);
    
/* Functions */
    int compute_ac(bool pacIN, int pac_strideIN, bool bpIN,
                    int bp_strideIN,bool spIN, int sp_strideIN,
                    bool wpIN,int wp_strideIN,std::string input_dir,
                    std::string output_dir,int num_samples=10,
                    int num_errors=10);
    
    int iPAC(int interval,std::string input_dir,std::string output_dir);
    
    int energy_conserve(double tol, int energy_stride, std::string input_dir,
                    std::string output_dir);

/* Mutators */
    
    /* All variables corresponding to those in set_system arguement are set
     after function call. is_sys_set = true after call */
    void set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                    double massIN, double betaIN, double beta_nuc_beadsIN,
                    double beta_elec_beadsIN, double alphaIN);
    
    /* All variables corresponding to those in set_time arguement are set
     after function call. is_time_set = true after call */
    void set_time(double dtIN, double total_timeIN);
    
    /* All variables corresponding to those in set_system arguement are set
     after function call. is_trajs_set = true after call and num_trajs_local
     is set as well*/
    void set_trajs(unsigned long long num_trajs_globalIN, std::string root_path);
    
private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    
    int nuc_beads, elec_beads, num_states; //system information
    double mass, beta, beta_nuc_beads, beta_elec_beads; //system information
    double alpha; //mapping variable pre-factor
    
    double dt, total_time; //time step size, total time [a.u.]
    unsigned long long num_trajs_global, num_trajs_local;
    
    //set to true once set_system, set_time, set_trajs, are called, otherwise false
    bool is_sys_set, is_time_set, is_trajs_set;
    
/* Functions */
    
    int pre_comp();

};

#endif

#ifndef auto_correlation_hpp
#define auto_correlation_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <string>
#include <fstream>
#include "mpi.h"

#include "M_Matrix.h"
#include "C_Matrix.h"
#include "dM_Matrix_dQ.hpp"
#include "theta_mixed.hpp"
#include "theta_mixed_dQ.hpp"
#include "theta_mixed_dElec.hpp"

#include "Forces_MTS.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "ABM_MVRPMD.hpp"

#include "aggregate.hpp"
#include "pop_estimators.hpp"
#include "trajs_io.hpp"

using namespace boost::numeric::ublas;

class auto_correlation{
   
public:
    
    auto_correlation(int my_id,int nun_procs, int root_proc);

    int compute(unsigned long long num_trajs_global,
                 unsigned long long num_trajs_local,
                 std::string input_dir,std::string output_dir,
                 int num_samples, int num_errors);

/*  Mutators */
    
    /* All variables corresponding to those in set_system arguement are set
     after function call. is_sys_set = true after call */
    void set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                    double massIN, double beta_nuc_beadsIN,
                    double beta_elec_beadsIN, double alphaIN);

    void request_calcs(bool pac, int pac_stride, bool bp, int bp_stride,
                       bool sp, int sp_stride, bool wp, int wp_stride);

    void set_time(double dtIN, double total_timeIN);
    
private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    
    double dt, total_time;
    
    int nuc_beads, elec_beads, num_states; //system information
    double mass, beta_nuc_beads, beta_elec_beads; //system information
    double alpha; //mapping variable pre-factor
    
    bool pac, bp, wp, sp;
    int pac_stride, bp_stride, sp_stride, wp_stride;
    
    bool is_sys_set, is_time_set, is_req_set;
    
/* Functions*/
    double compute_centroid(const vector<double> &Q);

    
};

#endif

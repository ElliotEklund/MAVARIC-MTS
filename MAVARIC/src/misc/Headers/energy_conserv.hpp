#ifndef energy_conserv_hpp
#define energy_conserv_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "GTerm.h"

#include "M_Matrix.h"
#include "C_Matrix.h"
#include "dM_Matrix_dQ.hpp"

#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_dQ.hpp"

#include "theta_mixed.hpp"
#include "theta_mixed_dQ.hpp"
#include "theta_mixed_dElec.hpp"

#include "Forces_MTS.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "mvrpmd_mixed_ham.hpp"
#include "ABM_MVRPMD.hpp"

#include <fstream>
#include <string>
#include "mpi.h"
#include "trajs_io.hpp"
#include <list>
#include <algorithm>

using namespace boost::numeric::ublas;

class energy_conserv{
    
public:
    
    energy_conserv(int my_id, int num_procs, int root_proc);
    
/* Functions */
    void compute(unsigned long long num_trajs_global,
                 unsigned long long num_trajs_local, double tol,
                 int energy_stride,std::string input_dir, std::string output_dir);

    
/* Mutators */
    void set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                    double massIN,double betaIN, double beta_nuc_beadsIN,
                    double beta_elec_beadsIN, double alphaIN);
    
    void set_time(double dtIN, double total_timeIN);

private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    int nuc_beads, elec_beads, num_states; //system data
    double mass, beta, beta_nuc_beads, beta_elec_beads; //system data
    double alpha; //mapping variable pre-factor
    bool is_sys_set, is_time_set;
    double dt, total_time; //time step, total running time
    
/* Functions */
    void write_broken(std::list<int> broken,unsigned long long num_trajs_local,
                      std::string file_root);
    
    void write_report(std::list <int> broken_trajectories,
                      unsigned long long num_trajs_global, double tol,
                      std::string output_dir);
};

#endif

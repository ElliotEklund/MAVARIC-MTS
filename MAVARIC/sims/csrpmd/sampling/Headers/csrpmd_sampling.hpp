#ifndef csrpmd_sampling_hpp
#define csrpmd_sampling_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "nr3.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"
#include "mpi.h"

#include <stdlib.h>
#include <time.h>
#include <string>

#include "csrpmd_forces.hpp"
#include "mv_forces_temp.hpp"
#include "RK4_MVRPMD.hpp"
#include "mpi_wrapper.hpp"

using namespace boost::numeric::ublas;

class csrpmd_sampling{
    
public:
    csrpmd_sampling(int my_id, int num_procs, int root_proc);
    
/* Functions */
    void run(std::string file_name);

/* Mutators */
    void set_sys_vars(int nuc_beads, int elec_beads, int num_states, double beta,
                      double mass);
    
    void set_sample_vars(int num_trajs,double decorr, double rate, double dt);
    
private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi information
    int nuc_beads, elec_beads, num_states; //system information
    double mass, beta; //mass, 1/kb T
    int num_trajs; //number of trajectories
    double dt, decorr, rate; //time step; decorrelation length; kick rate
    
    bool sys_set; //true if set_sys_vars has been called
    bool sample_set; //true if set_sample_vars has been called
    
/* Objects */
    Ran myRand; //NR3 random number generator
    
/* Functions */
    void init_Q(vector<double> &Q, double ss);
    
    void init_P(vector<double> &P);
    
    void init_elec(matrix<double> &x, matrix<double> &p);

    inline double step_dist(const double rn, double step_size);
    
    void save_trajs(vector<double> &v,std::string name);
    
    void save_trajs(matrix<double> &v,int size,int trajs, std::string name);

};

#endif

#ifndef sampling_mvrpmd_hpp
#define sampling_mvrpmd_hpp

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

#include "mpi_wrapper.hpp"
#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "SamplingHelper.h"
#include "system_step.hpp"
#include "elec_step.hpp"
#include "mvrpmd_mixed.hpp"
#include "mvrpmd_mixed_ham.hpp"

using namespace boost::numeric::ublas;

class sampling_mvrpmd{
    
public:
    sampling_mvrpmd(int my_id, int root_proc, int num_procs);
    
    /*
     Run main sampling routine
     nuc_ss, elec_ss: nuclear and electronic step sizes, respectively
     num_trajs: number of trajectories to be sampled
     decorr: decorrelation length between trajectories
     */
    int run(double nuc_ss, double x_ss, double p_xx,
             unsigned long long num_trajs,unsigned long long decorr);
    
    /*
     Set all variables in the argument; these are defined below under Model Data
     */
    void initialize_system(int nuc_beads_IN,int elec_beadsIN,int num_statesIN,
                           double massIN,double betaIN,double alphaIN);
    
    /* Set all variables in the argument; these are defined below under
     Sampling Data*/
    void initialize_files(bool readPSVIN,bool saveTrajsIN,
                          std::string rootFolderIN);
    
private:

/* Data */
    /* MPI Data*/
    int my_id, root_proc, num_procs;
    
    /* Model Data */
    int nuc_beads, elec_beads; //number of nuclear and electronic beads
    int num_states; //number of electronic states
    double beta; //1.0/kb T
    double mass; //system mass
    double alpha; //mapping variable prefactor
    
    /* Sampling Data */
    std::string rootFolder; //specifies folder to read and write data
    bool readPSV; //read in Phase Space Variables (PSV) if true
    bool saveTrajs;
    bool sys_set, files_set;
    
/* Objects */
    Ran myRand; //NR3 random number generator
    SamplingHelper helper;
    
/* Functions */
    /*
     Return a uniform random number between [-step_size,step_size]
     rn: random double
     step_size: sets range of uniform distribution
     */
    inline double step_dist(const double rn, double step_size);
    
    /* Save v to file name.
     v: vector to be saved
     name: name of file v is to be saved to*/
    void save_trajs(vector<double> &v,std::string name,int nuc_beads, int elec_beads,
                    int num_states,double betaIN,unsigned long long num_trajs_totalIN,
                    double decorrIN);
    
    /* Save v to file name.
     v: matrix to be saved
     size: number of states x number of beads in system
     num_trajs: number of trajectories v holds
     name: name of file v is to be saved to*/
    void save_trajs(matrix<double> &v,int size, unsigned long long num_trajs,
                    std::string name,int nuc_beads, int elec_beads, int num_states,
                    double betaIN,unsigned long long num_trajs_totalIN,double decorrIN);

    /* Initialize Q with random vaules.
     Q: vector of bead positions
     num_beads: number of beads in Q
     step_size: sets range of individual elements of Q; Q[i] ~ [-step_size,step_size]*/
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);
    
    /* Returns a uniform random number between [-step_size,step_size]
     rn: random double
     step_size: specifices range of uniform random distribution*/
    void gen_initElec(matrix<double> &v, int num_beads,
                      int num_states,double step_size);
    
    /* Return num_trajs_local, the correct number of trajectories for
     a given processor to run. */
    unsigned long long get_trajs_local(unsigned long long num_trajs_totalIN);

};

#endif

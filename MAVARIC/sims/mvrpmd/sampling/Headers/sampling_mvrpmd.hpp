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
    void run(double nuc_ss, double elec_ss,unsigned long long num_trajs,
                              unsigned long long decorr);
    
    /*
     Set all variables in the argument; these are defined below under Model Data
     */
    void initialize_system(int nuc_beads_IN,int elec_beadsIN,int num_statesIN,
                           double massIN,double betaIN);
    
    /* Set all variables in the argument; these are defined below under
     Sampling Data*/
    void initialize_files(bool readPSVIN, std::string rootFolderIN);

    
private:

/* Data */
    /* MPI Data*/
    int my_id, root_proc, num_procs;
    
    /* Model Data */
    int nuc_beads, elec_beads; //number of nuclear and electronic beads
    int num_states; //number of electronic states
    double beta; //1.0/kb T
    double mass; //system mass
    
    /* Sampling Data */
    std::string rootFolder; //specifies folder to read and write data
    bool readPSV; //read in Phase Space Variables (PSV) if true
    
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
    
    /*
     Return a uniform random number between [0,num_beads-1]
     rn: random Ullong; a nr3 data type
     num_beads: sets range of uniform distribution
     */
    inline int rand_bead(const Ullong rn, int num_beads);
    
    /*
     Return true if energy_prop passes Metropolis Hastings accpt/rejct criteria
     energy: currenty energy of system
     energy_prop: proposed energy of system
     */
    bool check_move(double energy, double energy_prop);
    
    void save_trajs(const vector<double> &v,std::string name);
    
    void save_trajs(const matrix<double> &v,int size, int num_trajs,
                    std::string name);

    
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);
    
    void gen_initElec(matrix<double> &v, int num_beads,
                      int num_states,double step_size);



    
};

#endif
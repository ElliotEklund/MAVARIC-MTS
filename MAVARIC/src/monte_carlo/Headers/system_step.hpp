#ifndef system_step_hpp
#define system_step_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "mvrpmd_mixed.hpp"
#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "mvrpmd_mixed_ham.hpp"

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class system_step{
    
public:
    system_step(int my_id,int num_procs, int root_proc,int num_beads,
                double betaIN);
    
    /* On average, try a move with every Q.
     energyIN: energy of system prior to calling step; will be changed by
               function call
     Q: bead positions of system prior to calling step; will be changed by function
     x: mapping variable prior to function call
     p: mapping variable prior to function call */
    void step(double energy,vector<double> &Q,const matrix<double> &x,
              const matrix<double> &p);
    
/* Mutators */
    void set_hamiltonian(mvrpmd_mixed &H_IN);
        
    void set_nuc_ss(double nuc_ssIN);
    
    void set_beta(double betaIN);
    
    void set_energy(double energyIN);
    
    double get_energy();

    unsigned long long get_steps_total();
    
    unsigned long long get_steps_accepted();
    
private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    Ran myRand; //NR3 random number generator
    int num_beads; //number of beads in system
    unsigned long long steps_total; //total number of steps taken
    unsigned long long steps_accepted; //total number of steps accepted
    double nuc_ss; //step size to take
    double beta; //1/kbT
    double energy; //energy of mc step

    vector<double> Q_prop; //proposed monte carlo move

/* Objects */
    //MVRPMD_MTS_Hamiltonian *H;
    mvrpmd_mixed *H;

/* Functions */
    
    /*
     Return a uniform random number between [0,num_beads-1]
     rn: random Ullong; a nr3 data type
     num_beads: sets range of uniform distribution
     */
    inline double step_dist(const double rn, double step_size);
    
    /*Return a uniform random number between [-step_size,step_size]
     rn: random double
     step_size: sets range of uniform distribution
     */
    inline int rand_bead(const Ullong rn, int num_beads);
    
    /* Return true is monte carlo criteria is met
     energy: energy before proposing move
     energy_prop: energy of proposed move*/
    bool check_move(double energy, double energy_prop);
};

#endif

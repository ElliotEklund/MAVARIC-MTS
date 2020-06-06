#ifndef elec_step_hpp
#define elec_step_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "mvrpmd_mixed.hpp"
#include "mvrpmd_mixed_ham.hpp"

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class elec_step{
    
public:
    elec_step(int my_id,int num_procs, int root_proc,int num_beads, int num_states,
              double betaIN);
    
    /* On average, try a move with every x mapping variable.
     energyIN: energy of system prior to calling step_x; will be changed by
               function call
     Q: bead positions of system prior to calling step_x
     x: mapping variable prior to function call; will be changed by function
     p: mapping variable prior to function call */
    void step_x(double energyIN, const vector<double> &Q, matrix<double> &x,
                const matrix<double> &p);
    
    /* On average, try a move with every p mapping variable.
     energyIN: energy of system prior to calling step_p; will be change by
               function call
     Q: bead positions of system prior to calling step_p
     x: mapping variable prior to function call
     p: mapping variable prior to function call; will be changed by function*/
    void step_p(double energyIN, const vector<double> &Q, const matrix<double> &x,
                matrix<double> &p);
    
/* Mutators */
    void set_hamiltonian(mvrpmd_mixed &H_IN);
        
    void set_beta(double betaIN);
    
    void set_energy(double energyIN);
    
    void set_ss(double x_ssIN, double p_ssIN);
    
    double get_energy();
    
    double get_x_steps_total();

    double get_x_steps_accpt();

    double get_p_steps_total();

    double get_p_steps_accpt();

private:
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    Ran myRand; //NR3 random number generator
    int num_beads; //number of beads in system
    int num_states; //number of electronic states
    double x_ss, p_ss; //x and p step sizes
    
    /* Number of x mapping variable steps tried and steps accepted*/
    unsigned long long x_steps_total, x_steps_accepted;
    
    /* Number of x mapping variable steps tried and steps accepted*/
    unsigned long long p_steps_total, p_steps_accepted;

    double beta; //1/kbT
    double energy; //energy of system before a move is proposed
    
    matrix<double> x_prop, p_prop; //propsed moves
    
/* Objects */
    mvrpmd_mixed *H; //Hamiltonian
    
/* Functions*/
    
    /*
     Return a uniform random number between [0,num_beads-1]
     rn: random Ullong; a nr3 data type
     num_beads: sets range of uniform distribution
     */
    inline int rand_bead(const Ullong rn, int num_beads);
    
    /*Return a uniform random number between [-step_size,step_size]
    rn: random double
    step_size: sets range of uniform distribution
    */
    inline double step_dist(const double rn, double step_size);
    
    /* Return true is monte carlo criteria is met
     energy: energy before proposing move
     energy_prop: energy of proposed move*/
    bool check_move(double energy, double energy_prop);
};

#endif

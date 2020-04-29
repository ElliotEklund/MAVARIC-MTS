#ifndef MonteCarloHelper_h
#define MonteCarloHelper_h

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "mpi.h"

using namespace boost::numeric::ublas;

class MonteCarloHelper{
    

public:
    
    MonteCarloHelper(std::string root, int my_id, int num_procs, int root_proc);
    
    // Print system Monte Carlo acceptance ratios after calling runMC.
    void print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt);
    
    //* Print electronic Monte Carlo acceptance ratios after calling runMC.  */
    void print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt);
    
    /* Print the MV-RPMD average energy after calling runMC from MonteCarlo */
    void print_avg_energy(double estimator_total, double sgn_total);
    
    /* Write estimator data to file "estimator" after calling runMC.  */
    void write_estimator(vector <double> estimator,int interval);
    
    void write_PSV(int nuc_beads,int elec_beads, int num_states,
                          vector<double> Q, matrix<double> x,matrix <double> p);
    
    /* Read in PSV.*/
    void read_PSV(int nuc_beads,int elec_beads, int num_states, vector<double> &Q, matrix<double> &x,
                         matrix<double> &p);
    
    void write_MC_data(double sgn_total,double estimator_total);
    
    void read_MC_data(double &sgn_totalGlobal, double &estimator_total);
    
    
private:
    std::string root;
    int my_id; //unique processor id
    int num_procs; //number of processors
    int root_proc; //root processor
    
};

#endif /* MonteCarloHelper_h */

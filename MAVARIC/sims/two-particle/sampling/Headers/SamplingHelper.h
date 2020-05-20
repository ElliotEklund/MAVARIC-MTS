#ifndef SamplingHelper_h
#define SamplingHelper_h

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

class SamplingHelper{
    
public:
    
    SamplingHelper(int my_id, int num_procs, int root_proc);
    
    void set_root(std::string rootIN);
    
    // Print system Monte Carlo acceptance ratios after calling runMC.
    void print_sys_accpt(unsigned long long sys_steps,
                         unsigned long long sys_steps_accpt, std::string sys_name);
    
    /* Read in PSV.*/
    void read_PSV(int nuc_beads,vector<double> &Q,std::string system);

private:
    std::string root;
    int my_id; //unique processor id
    int num_procs; //number of processors
    int root_proc; //root processor
    
};

#endif /* MonteCarloHelper_h */

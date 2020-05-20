#ifndef SYSTEM_H
#define SYSTEM_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

//#include "nr3.h"
//#include "ran.h"

#include <random>

using namespace boost::numeric::ublas;

class system{

public:
    
    system_move::system_move(int sys_beads, double sys_ss)
    
    vector<double> system_move(const vector<double> &Q);
    
private:
    
    vector<double> Q_prop;
    
    std::random_device myRand; //random number generator
    std::mt19937 mt; //specific implementation of random number generator
    std::uniform_real_distribution<double> sys_dist; //distribution of system steps
    std::uniform_int_distribution<int> rand_sys_bead; //generate random integer between 1 and sys_beads
}:

#endif

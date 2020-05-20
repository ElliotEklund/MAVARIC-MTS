#ifndef MOVE_SYSTEM_H
#define MOVE_SYSTEM_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <time.h>

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class move_system{

public:
    
    /* sys_beads: number of system ring polymer beads
       sys_ss: system step size for monte carlo move*/
    move_system(int sys_beads, double sys_ss);
    
    /*
     Return a vector of proposed moves
     Q: vector of system ring polymer positions of size
        sys_beads
     */
    vector<double> move(const vector<double> &Q);
    
private:
    
    vector<double> Q_prop; //contains proposed system moves
    double sys_ss; //system step size
    int sys_beads; //number of system ring polymer beads
    Ran myRand; //NR3 random number generator
    
    /* Given a random double, rn, return a random double between
     [-sys_ss, sys_ss].*/
    inline double sys_dist(const double rn);
    
    /* Given a random Ullong, rn, return a random system bead between
     [0:sys_beads)*/
    inline int rand_sys_bead(const Ullong rn);
    
};

#endif

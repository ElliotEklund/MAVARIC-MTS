#ifndef MOVE_BATH_H
#define MOVE_BATH_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <time.h>

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class move_bath{

public:
    
    /* sys_beads: number of system ring polymer beads
       sys_ss: system step size for monte carlo move*/
    move_bath(int bath_beads, int num_modes, vector<double> bath_ss);
    
    /*
     Return a vector of proposed moves
     Q: vector of system ring polymer positions of size
        sys_beads
     */
    matrix<double> move(const matrix<double> &Qbath);
    
private:
    
    matrix<double> Qbath_prop; //contains proposed system moves
    vector<double> bath_ss; //system step size
    int bath_beads; //number of system ring polymer beads
    int num_modes;
    Ran myRand; //NR3 random number generator
    
    /* Given a random double, rn, return a random double between
     [-sys_ss, sys_ss].*/
    inline double bath_dist(const double rn, const int mode);
    
    /* Given a random Ullong, rn, return a random system bead between
     [0:sys_beads)*/
    inline int rand_bath_bead(const Ullong rn);
    
    inline int rand_mode(const Ullong rn);
    
};

#endif

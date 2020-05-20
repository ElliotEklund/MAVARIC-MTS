#ifndef STATEINDEPPOT_H
#define STATEINDEPPOT_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <math.h>

using namespace boost::numeric::ublas;

class StateIndepPot{
    
public:
    StateIndepPot(int num_beads, double mass);
    
    /* Update ground_energy to reflect the state of Q
     Q is a vector of nuclear bead positions.*/
    void update_V0(const vector<double> &Q);
    
    /* Return ground_energy; update_ground_energy is called.
     Q is a vector of bead positions. */
    double& get_V0(const vector<double> &Q);
    
    /* Return ground_energy.
     Q is a vector of bead positions */
    double& get_V0();

private:
    
    /* Private data. */
    double mass; //system mass
    int num_beads;//number of nuclear beads
    double energy; //state independent energy
    
    vector<double> Q_SQ; //intermediate vector of bead positions squared
};

#endif

#include "move_bath.h"

move_bath::move_bath(int bath_beads, int num_modes, vector<double> bath_ss)
    :bath_beads(bath_beads),
     bath_ss(bath_ss),
     num_modes(num_modes),
     Qbath_prop(bath_beads,num_modes,0.0),
     myRand(time(NULL))
{}

matrix<double> move_bath::move(const matrix<double> &Qbath){
    
    int mode = rand_mode(myRand.int64()); //pick random mode
    
    for (int bead=0; bead<bath_beads; bead++) {
        //on average, move each bead
        double rand_bead = rand_bath_bead(myRand.int64());
        Qbath_prop(rand_bead,mode) = Qbath(rand_bead,mode) + bath_dist(myRand.int64(),mode);
    }
    
    return Qbath_prop;
}

inline double move_bath::bath_dist(const double rn, const int mode){
    return (rn * 2.0 * bath_ss(mode)) - bath_ss(mode);
}

inline int move_bath::rand_bath_bead(const Ullong rn){
    return rn % bath_beads;
}

inline int move_bath::rand_mode(const Ullong rn){
    return rn % num_modes;
}

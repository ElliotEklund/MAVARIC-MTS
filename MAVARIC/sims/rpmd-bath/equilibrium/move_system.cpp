#include "move_system.h"

move_system::move_system(int sys_beads, double sys_ss)
    :sys_beads(sys_beads),sys_ss(sys_ss),
     Q_prop(sys_beads,0.0),
     myRand(time(NULL))
{}

vector<double> move_system::move(const vector<double> &Q){
         
    double rand_bead = rand_sys_bead(myRand.int64());
    Q_prop = Q;
    Q_prop(rand_bead) = Q(rand_bead) + sys_dist(myRand.doub());
    return Q_prop;
}

inline double move_system::sys_dist(const double rn){
    return (rn * 2.0 * sys_ss) - sys_ss;
}

inline int move_system::rand_sys_bead(const Ullong rn){
    return rn % sys_beads;
}

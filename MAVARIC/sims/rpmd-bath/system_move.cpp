#include "system_move.h"

system_move::system_move(int sys_beads, double sys_ss)
    :Q_prop(sys_beads,0.0),
    mt(myRand()),
    nuc_dist(-sys_ss,sys_ss),
    rand_sys_bead(0,sys_beads - 1)


vector<double> system_move::step(const vector<double> &Q){
    int mcMove = rand_bead(mt);
    Q_prop(mcMove) = Q(mcMove) + nuc_dist(mt);
    return Q_prop;
}

#include "StateIndepPot.h"

StateIndepPot::StateIndepPot(int num_beads, double mass)
    :num_beads(num_beads), mass(mass),
     Q_SQ(num_beads)
{}
void StateIndepPot::update_V0(const vector<double> &Q){
    energy = 0.5 * inner_prod(Q,Q);
}
double& StateIndepPot::get_V0(const vector<double> &Q){
    
    update_V0(Q);
    return energy;
}
double & StateIndepPot::get_V0(){return energy;}

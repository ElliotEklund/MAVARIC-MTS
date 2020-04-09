#include "SpringEnergy.h"

SpringEnergy::SpringEnergy(int num_beads, double mass, double beta_num_beads)
    :num_beads(num_beads), mass(mass),
    beta_num_beads(beta_num_beads), Q_shift(num_beads),
    Q_diff(num_beads),
    preFactor(mass/(2.0*beta_num_beads*beta_num_beads))
{}

void SpringEnergy::update_springEnergy(const vector<double> &Q){
    
    update_Q_shift(Q);
    noalias(Q_diff) = Q - Q_shift;
    double prod = inner_prod(Q_diff, Q_diff);
    energy = preFactor*prod;
}

double& SpringEnergy::get_springEnergy(const vector<double> &Q){
    
    update_springEnergy(Q);
    return  energy;
}

double& SpringEnergy::get_springEnergy(){return  energy;}

void SpringEnergy::update_Q_shift(const vector<double> &Q){
    for (int bead=0; bead<num_beads; bead++) {
        Q_shift(bead) = Q((bead+1)%num_beads);
    }
}

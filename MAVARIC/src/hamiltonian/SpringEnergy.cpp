#include "SpringEnergy.h"

SpringEnergy::SpringEnergy(int num_beads, double mass, double beta_num_beads)
    :num_beads(num_beads), mass(mass),
    beta_num_beads(beta_num_beads),
    Q_diff(num_beads),
    preFactor(mass/(2.0*beta_num_beads*beta_num_beads)),
    W(num_beads,num_beads,2*num_beads)
{
    if (num_beads>1) {
        for (int i=0; i<num_beads; i++) {
            W(i,i) = 1.0;
            W(i,(i + 1)%num_beads) = - 1.0;
        }
    }
    else{
        W(0,0) = 0.0;
    }
}

void SpringEnergy::update_springEnergy(const vector<double> &Q){
    
    noalias(Q_diff) = prod(W,Q);
    double prod = inner_prod(Q_diff, Q_diff);
    energy = preFactor*prod;
}

double& SpringEnergy::get_springEnergy(const vector<double> &Q){
    
    update_springEnergy(Q);
    return  energy;
}

double& SpringEnergy::get_springEnergy(){return  energy;}

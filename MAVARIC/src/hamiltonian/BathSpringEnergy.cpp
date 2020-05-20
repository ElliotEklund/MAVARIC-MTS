#include "BathSpringEnergy.h"

BathSpringEnergy::BathSpringEnergy(int num_beads,int num_modes,double mass,
                                   double beta_num_beads)
    :num_modes(num_modes),
     prefactor(mass/(2.0*beta_num_beads*beta_num_beads)),
     W(num_beads,num_beads,2*num_beads),
     Q_diff(num_beads,num_modes,0),
     Q_diff_sq(num_beads,num_modes,0),
     ones_num_beads(num_beads,1.0),
     ones_num_modes(num_modes,1.0),
     temp_vec(num_modes,1.0)

{
    for (int i=0; i<num_beads; i++) {
        W(i,i) = 1.0;
        W(i,(i + 1)%num_beads) = - 1.0;
    }
}

void BathSpringEnergy::update_bathSpringEnergy(const matrix<double> &Q_bath){
    
    noalias(Q_diff) = prod(W,Q_bath);
    noalias(Q_diff_sq) = element_prod(Q_diff,Q_diff);
    noalias(temp_vec) = prod(ones_num_beads,Q_diff_sq);
    energy = prefactor*inner_prod(ones_num_modes,temp_vec);
}

double BathSpringEnergy::get_bathSpringEnergy(const matrix<double> &Q_bath){
    update_bathSpringEnergy(Q_bath);
    return energy;
}

double BathSpringEnergy::get_bathSpringEnergy(){return energy;}

#include "BathSpringEnergy.h"

BathSpringEnergy::BathSpringEnergy(int num_beads,int num_modes,double mass,
                                   double beta_num_beads)
    :num_modes(num_modes),V_spring(num_beads,mass,beta_num_beads)
{}

void BathSpringEnergy::update_bathSpringEnergy(const matrix<double> &Q_bath){
        
    double energy_temp = 0;
    
    for (int mode=0; mode<num_modes; mode++) {
        energy_temp += V_spring.get_springEnergy(column(Q_bath,mode));
    }
    energy = energy_temp;
}

double BathSpringEnergy::get_bathSpringEnergy(const matrix<double> &Q_bath){
    update_bathSpringEnergy(Q_bath);
    return energy;
}

double BathSpringEnergy::get_bathSpringEnergy(){return energy;}

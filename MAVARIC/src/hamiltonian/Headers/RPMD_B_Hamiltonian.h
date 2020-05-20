#ifndef RPMD_B_Hamiltonian_H
#define RPMD_B_Hamiltonian_H

#include "StateIndepPot.h"
#include "BathSpringEnergy.h"
#include "SpringEnergy.h"
#include "CouplingEnergy.h"

class RPMD_B_Hamiltonian{
    
public:

    RPMD_B_Hamiltonian(int sys_beads,double mass,double beta_sys_beads,
                       int bath_beads, int num_modes, double beta_bath_beads,
                       vector<double> cs, vector<double> ws);
    
    double get_energy(const vector<double> &Q,const matrix<double> &Qbath);

    double get_energy();
//    
//    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P,
//                          const matrix<double> &Qbath, const matrix<double> &Pbath);
    
private:
    int sys_beads;
    double mass;
    double beta_sys_beads;
    int bath_beads;
    int num_modes;
    double beta_bath_beads;
    
    StateIndepPot V0;
    SpringEnergy V_spring_sys;
    BathSpringEnergy V_spring_bath;
    CouplingEnergy V_couple;
    
    
};

#endif

#ifndef RPMD_Hamiltonian_hpp
#define RPMD_Hamiltonian_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"

class RPMD_Hamiltonian{
    
public:
    
    RPMD_Hamiltonian(double beta_num_beads,SpringEnergy &V_springIn,
                       StateIndepPot &V0_In);
    
    double get_energy(const vector<double> &Q);

    double get_energy();
    
    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P);
    
private:
    
    /* Private data*/
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    double one_beta_num_beads; // 1/beta_num_beads
    
};

#endif

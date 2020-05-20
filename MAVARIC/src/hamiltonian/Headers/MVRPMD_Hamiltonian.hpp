#ifndef MVRPMD_Hamiltonian_hpp
#define MVRPMD_Hamiltonian_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "Theta.h"
#include "GTerm.h"

class MVRPMD_Hamiltonian{
  
public:
    
    MVRPMD_Hamiltonian(double beta_num_beads,SpringEnergy &V_springIn, StateIndepPot &V0_In,
                       GTerm &GIn,Theta &thetaIn);
    
    /* Return energy corresponding to the state of Q, x, and p.
     This function implicitly updates all parts of the Hamiltonian.*/
    double get_energy(const vector<double> &Q, const matrix<double> &x,
                      const matrix<double> &p);
    
    /* Return energy. This function assumes that all parts of the
     Hamiltonian have already been updated.*/
    double get_energy();
    
private:
    
    /* Private data*/
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    GTerm * G;
    Theta * theta;
    double one_beta_num_beads; // 1/beta_num_beads
    
};
#endif

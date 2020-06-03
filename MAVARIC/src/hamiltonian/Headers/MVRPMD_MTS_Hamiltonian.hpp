#ifndef MVRPMD_MTS_Hamiltonian_hpp
#define MVRPMD_MTS_Hamiltonian_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "Theta_MTS.hpp"
#include "GTerm.h"
#include "mvrpmd_mixed.hpp"

class MVRPMD_MTS_Hamiltonian : public mvrpmd_mixed{

public:
    
    MVRPMD_MTS_Hamiltonian(double beta_num_beads,SpringEnergy &V_springIn,
                       StateIndepPot &V0_In,GTerm &GIn,Theta_MTS &thetaIn);

    double get_energy(const vector<double> &Q, const matrix<double> &x,
                      const matrix<double> &p);

    double get_energy();
    
    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P,
                          const matrix<double> &x,const matrix<double> &p);
    
    double get_sign();

private:

    /* Private data*/
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    GTerm * G;
    Theta_MTS * theta;
    double one_beta_num_beads; // 1/beta_num_beads

};

#endif

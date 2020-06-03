#ifndef mvrpmd_mixed_ham_hpp
#define mvrpmd_mixed_ham_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "theta_mixed.hpp"
#include "GTerm.h"
#include "mvrpmd_mixed.hpp"

class mvrpmd_mixed_ham : public mvrpmd_mixed{

public:
    
    mvrpmd_mixed_ham(double beta_num_beads,SpringEnergy &V_springIn,
                       StateIndepPot &V0_In,GTerm &GIn,theta_mixed &thetaIn);

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
    theta_mixed * theta;
    double one_beta_num_beads; // 1/beta_num_beads
};

#endif

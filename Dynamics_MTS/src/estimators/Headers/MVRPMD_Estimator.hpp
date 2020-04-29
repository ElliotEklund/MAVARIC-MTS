#ifndef MVRPMD_Estimator_hpp
#define MVRPMD_Estimator_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "Theta.h"
#include "dTheta_dBeta.hpp"


class MVRPMD_Estimator{
    
public:
    
    MVRPMD_Estimator(int num_beads, double beta_num_beads,
                     SpringEnergy &V_SpringIn, StateIndepPot &V0In,
                     Theta &thetaIn, dTheta_dBeta &dtheta_dBetaIn);
    
    /* Return energy estimator. This function assumes that all
     parts of the Hamiltonian have been updated. */
    double get_estimator();
    
private:
    
    /* Private data.*/
    double ONE_HALF_beta_num_beads; // 0.5 * (1/beta_num_beads)
    double ONE_num_beads; //1.0/num_beads
    
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    Theta * theta;
    dTheta_dBeta * dtheta_dBeta;
    
};

#endif

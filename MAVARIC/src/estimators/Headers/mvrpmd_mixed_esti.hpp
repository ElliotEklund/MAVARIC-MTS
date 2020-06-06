#ifndef mvrpmd_mixed_esti_hpp
#define mvrpmd_mixed_esti_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
//#include "Theta_MTS.hpp"
//#include "dTheta_MTS_dBeta.hpp"

#include "theta_mixed.hpp"
#include "theta_mixed_dBeta.hpp"

class mvrpmd_mixed_esti{
    
public:
    
    mvrpmd_mixed_esti(int num_beads, double beta_num_beads,
                      SpringEnergy &V_SpringIn, StateIndepPot &V0In,
                      theta_mixed &thetaIN, theta_mixed_dBeta &dthetaIN);
    
    double get_estimator();
    
    double get_estimator(const vector<double> &Q,const matrix<double> &x,
                         const matrix<double> &p);
    
private:
    
/* Data */
    double ONE_HALF_beta_num_beads; // 0.5 * (1/beta_num_beads)
    double ONE_num_beads;
    
/* Objects*/
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    theta_mixed * theta;
    theta_mixed_dBeta * dtheta;
};

#endif

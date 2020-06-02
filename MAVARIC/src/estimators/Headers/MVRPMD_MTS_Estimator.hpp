#ifndef MVRPMD_MTS_Estimator_hpp
#define MVRPMD_MTS_Estimator_hpp


#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "Theta_MTS.hpp"
#include "dTheta_MTS_dBeta.hpp"


class MVRPMD_MTS_Estimator{
    
public:
    
    MVRPMD_MTS_Estimator(int num_beads, double beta_num_beads,
                     SpringEnergy &V_SpringIn, StateIndepPot &V0In,
                     Theta_MTS &thetaMTSIn, dTheta_MTS_dBeta &dthetaMTS_dBetaIn);
    
    double get_estimator();
    
    double get_estimator(const vector<double> &Q,const matrix<double> &x,
                         const matrix<double> &p);
    
private:
    
    double ONE_HALF_beta_num_beads; // 0.5 * (1/beta_num_beads)
    double ONE_num_beads;
    
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    Theta_MTS * thetaMTS;
    dTheta_MTS_dBeta * dthetaMTS_dBeta;
    
};


#endif

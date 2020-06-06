#include "MVRPMD_MTS_Estimator.hpp"

MVRPMD_MTS_Estimator::MVRPMD_MTS_Estimator(int num_beads, double beta_num_beads,
                                           SpringEnergy &V_SpringIn,
                                           StateIndepPot &V0In,
                                           Theta_MTS &thetaMTSIn,
                                           dTheta_MTS_dBeta &dthetaMTS_dBetaIn)
    :ONE_HALF_beta_num_beads(0.5/beta_num_beads),
     ONE_num_beads(1.0/num_beads)
{
    V_spring = &V_SpringIn;
    V0 = &V0In;
    thetaMTS = &thetaMTSIn;
    dthetaMTS_dBeta = &dthetaMTS_dBetaIn;
}

double MVRPMD_MTS_Estimator::get_estimator(){
    
    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    double thetaMTS_temp = thetaMTS->get_theta();
    double sgnThetaMTS_temp = thetaMTS->get_signTheta();
    double dthetaMTS_dBeta_temp = dthetaMTS_dBeta->get_dTheta_MTS_dBeta();
    
    return (ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp - V_spring_temp) -
            dthetaMTS_dBeta_temp/thetaMTS_temp)*sgnThetaMTS_temp;
}

double MVRPMD_MTS_Estimator::get_estimator(const vector<double> &Q,
                                           const matrix<double> &x,
                                           const matrix<double> &p){
    
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    double thetaMTS_temp = thetaMTS->get_theta(Q,x,p);
    double sgnThetaMTS_temp = thetaMTS->get_signTheta();
    double dthetaMTS_dBeta_temp = dthetaMTS_dBeta->get_dTheta_MTS_dBeta();
     
    return (ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp - V_spring_temp)
            - dthetaMTS_dBeta_temp/thetaMTS_temp)*sgnThetaMTS_temp;
}

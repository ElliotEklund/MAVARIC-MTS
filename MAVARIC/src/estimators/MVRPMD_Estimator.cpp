#include "MVRPMD_Estimator.hpp"

MVRPMD_Estimator::MVRPMD_Estimator(int num_beads, double beta_num_beads,                                            SpringEnergy &V_SpringIn, StateIndepPot &V0In,
                                    Theta &thetaIn, dTheta_dBeta &dtheta_dBetaIn)
    :ONE_HALF_beta_num_beads(0.5/beta_num_beads),
     ONE_num_beads(1.0/num_beads)
{
    V_spring = &V_SpringIn;
    V0 = &V0In;
    theta = &thetaIn;
    dtheta_dBeta = &dtheta_dBetaIn;
}

double MVRPMD_Estimator::get_estimator(){
    
    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    double theta_temp = theta->get_theta();
    double sgnTheta_temp = theta->get_signTheta();
    double dtheta_dBeta_temp = dtheta_dBeta->get_dTheta_dBeta();

    return (ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp - V_spring_temp) - dtheta_dBeta_temp/theta_temp)*sgnTheta_temp;
}

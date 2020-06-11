#include "mvrpmd_mixed_ham.hpp"

mvrpmd_mixed_ham::mvrpmd_mixed_ham(double beta_num_beads,SpringEnergy
                                       &V_springIn,StateIndepPot &V0_In,
                                       GTerm &GIn, theta_mixed &thetaIn)
    :one_beta_num_beads(1.0/beta_num_beads)
{
    V_spring = &V_springIn;
    V0 = &V0_In;
    G = &GIn;
    theta = &thetaIn;
}
double mvrpmd_mixed_ham::get_energy(const vector<double> &Q,const matrix<double> &x,
                                      const matrix<double> &p){
    
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    double G_temp = G->get_gTerm(x, p);
    double theta_temp = theta->get_theta(Q, x, p);
        
    return V_spring_temp + V0_temp + one_beta_num_beads*(G_temp - log(fabs(theta_temp)));
}

double mvrpmd_mixed_ham::get_energy(){
    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    double G_temp = G->get_gTerm();
    double theta_temp = theta->get_theta();
        
    return V_spring_temp + V0_temp + one_beta_num_beads*(G_temp - log(fabs(theta_temp)));
}

double mvrpmd_mixed_ham::get_energy_dyn(double mass,
                                              const vector<double> &Q, const vector<double> &P,
                                              const matrix<double> &x,const matrix<double> &p){
    
    double kin_energy = 0.5 * inner_prod(P,P)/mass;
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    double G_temp = G->get_gTerm(x, p);
    double theta_temp = theta->get_theta();
    
    return kin_energy + V_spring_temp + V0_temp +
           one_beta_num_beads*(G_temp - log(fabs(theta_temp)));
}
double  mvrpmd_mixed_ham::get_sign(){
    return theta->get_signTheta();
}

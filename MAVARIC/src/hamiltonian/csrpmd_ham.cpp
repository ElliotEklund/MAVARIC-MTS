#include "csrpmd_ham.hpp"

csrpmd_ham::csrpmd_ham(int nuc_beads, int elec_beads,double beta_num_beads,
                       SpringEnergy &V_springIn,StateIndepPot &V0_In,
                       sc_potential &Vsc_IN)
    :nuc_beads(nuc_beads),
     elec_beads(elec_beads),
     one_beta_num_beads(1.0/beta_num_beads)
{
    V_spring = &V_springIn;
    V0 = &V0_In;
    Vsc = &Vsc_IN;
}
double csrpmd_ham::get_energy(const vector<double> &Q,const matrix<double> &x,
                          const matrix<double> &p){
    
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    double Vsc_temp = Vsc->get_Vsc(Q,x,p);
    return V_spring_temp + V0_temp + Vsc_temp;
}
double csrpmd_ham::get_energy(){
 
    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    double Vsc_temp = Vsc->get_Vsc();
    return V_spring_temp + V0_temp + Vsc_temp;
}
double csrpmd_ham::get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P,
                                  const matrix<double> &x,const matrix<double> &p){
    
    double kin_energy = 0.5 * inner_prod(P,P)/mass;
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    double Vsc_temp = Vsc->get_Vsc(Q,x,p);
    return kin_energy + V_spring_temp + V0_temp + Vsc_temp;
}

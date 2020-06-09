#include "RPMD_Hamiltonian.hpp"

RPMD_Hamiltonian::RPMD_Hamiltonian(double beta_num_beads,SpringEnergy &V_springIn,
                                   StateIndepPot &V0_In)
    :one_beta_num_beads(1.0/beta_num_beads)
{
    V_spring = &V_springIn;
    V0 = &V0_In;
}
double RPMD_Hamiltonian::get_energy(const vector<double> &Q){
    
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    
    return V_spring_temp + V0_temp;
}
double RPMD_Hamiltonian::get_energy(){
 
    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    return V_spring_temp + V0_temp;
}
double RPMD_Hamiltonian::get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P){
    
    double kin_energy = 0.5 * inner_prod(P,P)/mass;
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    
    return kin_energy + V_spring_temp + V0_temp;
}

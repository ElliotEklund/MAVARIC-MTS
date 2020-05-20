#include "RPMD_B_Hamiltonian.h"

RPMD_B_Hamiltonian::RPMD_B_Hamiltonian(int sys_beads,double mass, double beta_sys_beads,
                                       int bath_beads, int num_modes, double beta_bath_beads,
                                       vector<double> cs, vector<double> ws)
    :sys_beads(sys_beads),mass(mass),beta_sys_beads(beta_sys_beads),
     bath_beads(bath_beads), num_modes(num_modes), beta_bath_beads(beta_bath_beads),
     V0(sys_beads,mass),
     V_spring_sys(sys_beads,mass,beta_sys_beads),
     V_spring_bath(bath_beads,num_modes,mass,beta_bath_beads),
     V_couple(num_modes,sys_beads,bath_beads,mass, cs, ws)
{
    int i = 0;
}

double RPMD_B_Hamiltonian::get_energy(const vector<double> &Q,const matrix<double> &Qbath){

    double E_spring_sys = V_spring_sys.get_springEnergy(Q);
    double E_spring_bath = V_spring_bath.get_bathSpringEnergy(Qbath);
    double E_v0 = V0.get_V0(Q);
    double E_couple = V_couple.get_couplingEnergy(Q,Qbath);
    
    return E_spring_sys + E_spring_bath + E_v0 + E_couple;
}

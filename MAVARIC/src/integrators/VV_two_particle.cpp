#include "VV_two_particle.hpp"

VV_two_particle::VV_two_particle(int nuc_beads1, int nuc_beads2,
                                 double mass, double beta, double dt)
    :dt(dt), HALF_dt(0.5*dt),
     F(nuc_beads1,nuc_beads2,mass,beta),
     P1_half(nuc_beads1,0),P2_half(nuc_beads1,0)
{}
void VV_two_particle::step(vector<double> &Q1,vector<double> &Q2,
                           vector<double> &P1,vector<double> &P2){
    
    F.update_dHdQ1(Q1,Q2);
    F.update_dHdQ2(Q1,Q2);
    
    noalias(P1_half) = P1 - HALF_dt*F.get_dHdQ1();
    noalias(P2_half) = P2 - HALF_dt*F.get_dHdQ2();
        
    F.update_dHdP1(P1_half);
    F.update_dHdP2(P2_half);

    Q1 = Q1 + dt*F.get_dHdP1();
    Q2 = Q2 + dt*F.get_dHdP2();
        
    F.update_dHdQ1(Q1,Q2);
    F.update_dHdQ2(Q1,Q2);
    
    noalias(P1) = P1_half - HALF_dt*F.get_dHdQ1();
    noalias(P2) = P2_half - HALF_dt*F.get_dHdQ2();

}

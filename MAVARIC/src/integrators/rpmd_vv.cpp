#include "rpmd_vv.hpp"

rpmd_vv::rpmd_vv(int nuc_beads,double mass, double beta, double dt)
    :dt(dt), HALF_dt(0.5*dt),
     F(nuc_beads,mass,beta),
     P_half(nuc_beads,0)
{}
void rpmd_vv::step(vector<double> &Q,vector<double> &P){
    
    F.update_dHdQ(Q);
    noalias(P_half) = P - HALF_dt*F.get_dHdQ();
        
    F.update_dHdP(P_half);
    Q = Q + dt*F.get_dHdP();
        
    F.update_dHdQ(Q);
    noalias(P) = P_half - HALF_dt*F.get_dHdQ();
}

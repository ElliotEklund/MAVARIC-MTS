#include "mvrpmd_special.hpp"

mvrpmd_special::mvrpmd_special(Forces_MTS *F_In,int nuc_beads,
                               int elec_beads,int num_states,double dt)
    :dt(dt),
     P_half(nuc_beads,0),
     x_next(elec_beads,num_states,0),
     p_next(elec_beads,num_states,0),
     x_half(elec_beads,num_states,0),
     p_half(elec_beads,num_states,0)
     
{
    F = F_In;
    dt_half = 0.5*dt;
    dt_tenth = 0.1*dt;
}
void mvrpmd_special::step(vector<double> &Q, vector<double> &P,
          matrix<double> &x, matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    P_half = P - dt_half * F->get_dHdQ();

    for (int i=0; i<5; i++) {
        x_next = x + dt_tenth * F->get_dHdp();
        p_next = p - dt_tenth * F->get_dHdx();
        F->update_Forces(Q,P,x_next,p_next);
        x = x_next;
        p = p_next;
    }

    x_half = x;
    p_half = p;

    for (int i=0; i<5; i++) {
        x_next = x + dt_tenth * F->get_dHdp();
        p_next = p - dt_tenth * F->get_dHdx();
        F->update_Forces(Q,P,x_next,p_next);
        x = x_next;
        p = p_next;
    }

    F->update_Forces(Q,P_half,x,p);
    Q = Q + dt*F->get_dHdP();

    F->update_Forces(Q,P_half,x_half,p_half);
    P = P_half - dt_half*F->get_dHdQ();
}

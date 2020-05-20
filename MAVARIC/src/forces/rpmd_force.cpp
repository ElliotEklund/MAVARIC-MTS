#include "rpmd_force.hpp"

rpmd_force::rpmd_force(int nuc_beads,double mass, double beta)
    :nuc_beads(nuc_beads),
     dHdQ(nuc_beads,0.0),
     dHdP(nuc_beads,0.0),
     dVspring_dQ(nuc_beads,mass,beta/nuc_beads),
     dVspring_dQ_vec(nuc_beads,0),
     dV0_dQ_vec(nuc_beads,0)

{
    c = 1.0/mass;
}

void rpmd_force::update_dHdP(const vector<double> &P){
    noalias(dHdP) = nuc_beads * c * P;
}
void rpmd_force::update_dHdQ(const vector<double> &Q){
    noalias(dVspring_dQ_vec) = dVspring_dQ.get_dSpring_dQ2(Q);
    noalias(dV0_dQ_vec) = dV0_dQ(Q);
    noalias(dHdQ) = (dVspring_dQ_vec + dV0_dQ_vec)/nuc_beads;
}
vector<double> rpmd_force::dV0_dQ(const vector<double> &Q){
    vector<double> QQ = element_prod(Q,Q);
    vector<double> QQQ = element_prod(QQ,Q);
    return Q + 0.3*QQ + 0.04*QQQ;
}

const vector<double> & rpmd_force::get_dHdQ(){return dHdQ;}

const vector<double> & rpmd_force::get_dHdP(){return dHdP;}

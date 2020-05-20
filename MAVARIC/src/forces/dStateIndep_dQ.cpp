#include "dStateIndep_dQ.hpp"

dStateIndep_dQ::dStateIndep_dQ(int nuc_beads,double mass)
    :nuc_beads(nuc_beads), mass(mass),

     force(nuc_beads,0.0)
{}

const vector<double> & dStateIndep_dQ::get_dStateIndep_dQ(const vector<double> &Q){
    force = Q;
    return force;
}

#include "dSpring_dQ.hpp"

dSpring_dQ::dSpring_dQ(int nuc_beads,double mass, double beta_nuc_beads)
    :nuc_beads(nuc_beads),mass(mass),beta_nuc_beads(beta_nuc_beads),
     coeff(mass/(beta_nuc_beads*beta_nuc_beads)),
     W(nuc_beads,nuc_beads),
     force(nuc_beads,0.0)
{
    for (int i=0; i<nuc_beads; i++) {
        W(i,i) = 2.0;
        W(i,(i+1)%nuc_beads) = -1.0;
        W(i,(i+nuc_beads-1)%nuc_beads) = -1.0;
    }
}

/* Return the derivative of the spring term w.r.t to Q*/
const vector<double> & dSpring_dQ::get_dSpring_dQ(const vector<double> &Q){
    
    int i_up = 0; //index one bead up
    int i_down = 0; //index one bead down
    
    for(int bead=0; bead<nuc_beads; bead++){
        
        i_up = (bead + 1) % nuc_beads;
        i_down = (bead + nuc_beads - 1) % nuc_beads;
        
        force(bead) = 2.0*Q(bead) - Q(i_up) - Q(i_down);
    }

    force = coeff * force;
    return force;
}

/* Return the derivative of the spring term w.r.t to Q*/
const vector<double> & dSpring_dQ::get_dSpring_dQ2(const vector<double> &Q){
    
    noalias(force) = coeff*prod(W,Q);
    return force;
}

#include "CouplingEnergy.h"

CouplingEnergy::CouplingEnergy(int bath_modes, int sys_beads, int bath_beads,
                               double mass, vector<double> cs, vector<double> ws)
    :bath_modes(bath_modes), sys_beads(sys_beads), bath_beads(bath_beads),
     r(sys_beads/bath_beads),mass(mass),
     cs(cs), ws(ws),
     wwm(bath_modes), c_wwm(bath_modes),
     W(bath_beads,sys_beads),
     ones_bath_beads(bath_beads,1.0),
     Q_sub(bath_beads,0.0),
     cQ_sub(bath_beads,bath_modes),
     dif(bath_beads,bath_modes,0.0),
     dif_sq(bath_beads,bath_modes,0.0),
     sum_dif_sq(bath_modes,0.0)

{
    for (int mode=0; mode<bath_modes; mode++) {
        wwm(mode) = 0.5*ws[mode]*ws[mode]*mass;
        c_wwm(mode) = cs(mode)/(wwm(mode));
    }
    
    for (int i=0; i<bath_beads; i++) {
        int s = i*r;
        for (int j=0; j<r; j++) {
            W(i,i+j) = 1.0;
        }
    }
}

void CouplingEnergy::update_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath){
    
    noalias(Q_sub) = prod(W,Q);
    noalias(cQ_sub) = outer_prod(Q_sub,c_wwm);
    noalias(dif) = Q_bath - cQ_sub;
    noalias(dif_sq) = element_prod(cQ_sub,cQ_sub);
    noalias(sum_dif_sq) = prod(ones_bath_beads,dif_sq);
    energy = inner_prod(wwm,sum_dif_sq);
}

double CouplingEnergy::get_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath){
    update_couplingEnergy(Q, Q_bath);
    return energy;
}

double CouplingEnergy::get_couplingEnergy(){return energy;}

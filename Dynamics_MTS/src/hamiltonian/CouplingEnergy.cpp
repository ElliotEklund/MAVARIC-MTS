#include "CouplingEnergy.h"

CouplingEnergy::CouplingEnergy(int bath_modes, int num_beads, double mass,
                               vector<double> cs, vector<double> ws)
    :bath_modes(bath_modes), num_beads(num_beads), mass(mass),
     cs(cs), ws(ws), wwm(bath_modes), c_wwm(bath_modes),
     diff_sq(num_beads)
{
    for (int mode=0; mode<bath_modes; mode++) {
        wwm(mode) = 0.5*ws[mode]*ws[mode]*mass;
        c_wwm(mode) = cs(mode)/(wwm(mode));
    }
}

void CouplingEnergy::update_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath){
 
    energy = 0;

    for (int mode=0; mode<bath_modes; mode++) {
        diff_sq = column(Q_bath,mode) - c_wwm(mode)*Q;
        std::transform(diff_sq.begin(), diff_sq.end(), diff_sq.begin(), square());
        energy += wwm(mode) * sum(diff_sq);
    }
}

double CouplingEnergy::get_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath){
    
    update_couplingEnergy(Q, Q_bath);
    return energy;
}

double CouplingEnergy::get_couplingEnergy(){return energy;}

double CouplingEnergy::square::operator() (double x) const { return pow(x, 2); }


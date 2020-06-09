#ifndef csrpmd_ham_hpp
#define csrpmd_ham_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "sc_potential.hpp"

class csrpmd_ham{
    
public:
    
    csrpmd_ham(int nuc_beads, int elec_beads, double beta_num_beads,
               SpringEnergy &V_springIn,StateIndepPot &V0_In, sc_potential & Vsc_IN);

/* Mutators */
    double get_energy(const vector<double> &Q,const matrix<double> &x,
                      const matrix<double> &p);

    double get_energy();
    
    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P,
                         const matrix<double> &x, const matrix<double> &p);
    
private:
    
/* Data*/
    int nuc_beads, elec_beads; //nuclear and electronic ring polymer beads
    double one_beta_num_beads; // 1/beta_num_beads
    
/* Objects */
    SpringEnergy * V_spring;
    StateIndepPot * V0;
    sc_potential * Vsc;
};

#endif

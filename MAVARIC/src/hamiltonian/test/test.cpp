#include <iostream>
#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "rpmd_vv.hpp"
#include "RK4_MVRPMD.hpp"
#include "Forces_MTS.hpp"
#include "mvrpmd_special.hpp"
#include "sc_potential.hpp"
#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "csrpmd_ham.hpp"
#include "csrpmd_forces.hpp"

using namespace boost::numeric::ublas;

int main(int argc, char ** argv){

    int nuc_beads, elec_beads, num_states;
    nuc_beads = 1;
    elec_beads = 1;
    num_states = 2;
    double beta = 1.0;
    double mass = 1.0;
    double dt = 0.1;
    
    vector<double> Q(nuc_beads,0);
    vector<double> P(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0);
    matrix<double> p(elec_beads,num_states,0);

    for(int bead=0; bead<nuc_beads; bead++){
        Q(bead) = bead - 3;
        P(bead) = bead*2.78;
    }

    for(int bead=0;bead<elec_beads;bead++){
        for(int state=0;state<num_states;state++){
            x(bead,state) = bead - 10*state;
            p(bead,state) = bead + state;
        }
    }

    sc_potential Ve(nuc_beads,elec_beads,num_states);
    SpringEnergy Vspring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);

    csrpmd_ham H(nuc_beads,elec_beads,beta/nuc_beads,Vspring,V0,Ve);
    double energy = H.get_energy(Q,x,p);

    csrpmd_forces F(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads);
    //mv_forces_temp *F_temp = &F;
    
    F.update_Forces(Q,P,x,p);
    F.print_dHdQ();
    F.print_dHdx();
    F.print_dHdp();

    //RK4_MVRPMD rk4(F_temp,nuc_beads,elec_beads,num_states,dt);
 
    //rk4.take_step(Q,P,x,p);

    return 0;
}

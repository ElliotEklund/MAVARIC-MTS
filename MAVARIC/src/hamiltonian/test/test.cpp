#include <iostream>
#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "mvrpmd_mixed_ham.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "theta_mixed_dQ.hpp"


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
            x(bead,state) = 0.1*bead - 0.2*state;
            p(bead,state) = 0.1*bead + 0.3*state;
        }
    }
    
    SpringEnergy Vspring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    C_Matrix C(elec_beads,num_states);
    M_Matrix M(num_states,elec_beads,beta/elec_beads);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
    mvrpmd_mixed_ham H(beta/nuc_beads,Vspring,V0,G,theta);

    H.get_energy(Q,x,p);
    
    dM_Matrix_dQ M_dQ(elec_beads,num_states,beta/elec_beads,M);
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,M_dQ);
    theta_mixed_dElec theta_dElec(num_states,elec_beads,C,M);
    
    mvrpmd_mixed_forces F(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
                          theta,theta_dQ,theta_dElec);
    
    F.update_Forces(Q,P,x,p);

    return 0;
}

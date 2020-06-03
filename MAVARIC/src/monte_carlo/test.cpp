#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "system_step.hpp"
#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "Theta_MTS.hpp"
#include "theta_mixed.hpp"
#include "GTerm.h"
#include "mvrpmd_mixed.hpp"
#include "mvrpmd_mixed_ham.hpp"

using namespace boost::numeric::ublas;

int main(){

    std::cout << "Hello, world!" << std::endl;

    int my_id = 0;
    int num_procs = 1;
    int root_proc = 0;
    int nuc_beads = 6;
    double beta = 1.0;
    double mass = 1.0;
    int num_states = 2;
    double elec_beads = 3;

    double ss = 0.5;
    double energy = 1000.0;

    system_step stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    stepper.set_nuc_ss(ss);


        /* Assemble Hamiltonian and Estimator*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    C_Matrix C(elec_beads,num_states);
    M_Matrix M(num_states,nuc_beads,beta/elec_beads);
    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M);
    Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);
    theta_mixed theta_2(num_states,nuc_beads,elec_beads,C,M);

    MVRPMD_MTS_Hamiltonian H(beta/nuc_beads,V_spring,V0,G,thetaMTS);
    mvrpmd_mixed_ham H_2(beta/nuc_beads,V_spring,V0,G,theta_2);
    
    stepper.set_hamiltonian(H_2);

    vector<double> Q(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0),p(elec_beads,num_states,0);

    for (int i=0; i<nuc_beads; i++){
        Q(i) = i - 3.0;
    }

    for(int i=0; i<elec_beads; i++){
        for(int j=0; j<num_states; j++){
            x(i,j) = i+j;
            p(i,j) = i-j;
        }
    }


    stepper.step(energy,Q,x,p);
    std::cout << stepper.get_energy() << std::endl;
    return 0;
}

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
#include "C_Matrix.h"

#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_MTS_dQ.hpp"
#include "dM_Matrix_dQ.hpp"

#include "dTheta_MTS_dQ.hpp"
#include "dTheta_MTS_dElec.hpp"
#include "theta_mixed_dQ.hpp"

#include "MVRPMD_MTS_Estimator.hpp"
#include "mvrpmd_mixed_esti.hpp"

#include "mvrpmd_mixed_forces.hpp"
#include "Forces_MTS.hpp"
#include "theta_mixed_dElec.hpp"

#include "RK4_MVRPMD.hpp"
#include "ABM_MVRPMD.hpp"

using namespace boost::numeric::ublas;

int main(){

    int my_id = 0;
    int num_procs = 1;
    int root_proc = 0;
    int nuc_beads = 4;
    double beta = 1.0;
    double mass = 1.0;
    int num_states = 2;
    double elec_beads = 4;

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
    dM_Matrix_dQ dMdQ (nuc_beads, num_states, beta/elec_beads, M);
    
    theta_mixed theta_2(num_states,nuc_beads,elec_beads,C,M);
    
    theta_mixed_dBeta dtheta_2(elec_beads,num_states,beta/elec_beads,C,M);
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ);
    theta_mixed_dElec theta_dElec(num_states, elec_beads, C, M);
    
    mvrpmd_mixed_ham H_2(beta/nuc_beads,V_spring,V0,G,theta_2);
    mvrpmd_mixed_esti Esti_2(nuc_beads,beta/nuc_beads,V_spring,V0,theta_2,
                                dtheta_2);

    
    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M);
    dM_Matrix_MTS_dQ dM_MTS_dQ(nuc_beads, elec_beads, num_states, dMdQ);
    dTheta_MTS_dQ dThetadQ(num_states, nuc_beads, elec_beads, C, M_MTS, dM_MTS_dQ);
    dTheta_MTS_dElec dThetadElec(num_states, elec_beads, C, M_MTS);
    Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);

    //MVRPMD_MTS_Hamiltonian H(beta/nuc_beads,V_spring,V0,G,thetaMTS);
    
    dTheta_MTS_dBeta dtheta_1(nuc_beads,elec_beads,num_states,beta/elec_beads,C,M,M_MTS);
    
//    MVRPMD_MTS_Estimator Esti_1(nuc_beads,beta/nuc_beads,V_spring,V0,thetaMTS,
//                                dtheta_1);
    
    mvrpmd_mixed_forces F(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
                          theta_2,theta_dQ,theta_dElec);
    
    Forces_MTS F2(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
                            thetaMTS,dThetadQ,dThetadElec);
    

    double dt = 0.01;
    ABM_MVRPMD myStepper(F,dt,num_states,nuc_beads,elec_beads);
    ABM_MVRPMD myStepper2(F2,dt,num_states,nuc_beads,elec_beads);


    vector<double> Q(nuc_beads,0);
    vector<double> P(nuc_beads,0);

    matrix<double> x(elec_beads,num_states,0),p(elec_beads,num_states,0);

    for (int i=0; i<nuc_beads; i++){
        Q(i) = i - 3.0;
        P(i) = i + 1.0;
    }

    for(int i=0; i<elec_beads; i++){
        for(int j=0; j<num_states; j++){
            x(i,j) = i+j;
            p(i,j) = i-j;
        }
    }
    
    std::cout << Q << std::endl;
    std::cout << P << std::endl;
    std::cout << x << std::endl;
    std::cout << p << std::endl;
    
    myStepper2.initialize_rk4(Q,P,x,p);
    myStepper2.take_step(Q,P,x,p);
    
    std::cout << Q << std::endl;
    std::cout << P << std::endl;
    std::cout << x << std::endl;
    std::cout << p << std::endl;



//    C.update_C_vec(x,p);
//    M.update_M_vec(Q);
//    dMdQ.update_dM_dQ_vec(Q);
//    theta_dQ.update_theta_dQ(Q);
//
//    M_MTS.update_M_MTS_vec(Q);
//    dM_MTS_dQ.update_dM_MTS_dQ_vec(Q);
//    dThetadQ.update_dTheta_MTS_dQ_vec(Q);
//
//    std::cout << theta_dQ.get_theta_dQ_vec() << std::endl;
//    std::cout << dThetadQ.get_dThetaMTS_dQ_vec() << std::endl;
    
    return 0;
}

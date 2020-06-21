#include <iostream>
#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "StateIndepPot.h"
#include "theta_mixed_dQ.hpp"
#include "mvrpmd_mixed_forces.hpp"

#include "Theta_MTS.hpp"
#include "dTheta_MTS_dQ.hpp"
#include "dTheta_MTS_dElec.hpp"
#include "Forces_MTS.hpp"

#include "MVRPMD_MTS_Hamiltonian.hpp"
#include "mvrpmd_mixed_ham.hpp"
#include "ABM_MVRPMD.hpp"
#include <limits>
#include "functions.hpp"

vector<double> get_Q(){
    
    vector<double> Q(6,0);
    
    Q(0) = 0.065518;
    Q(1) =0.471757;
    Q(2) =0.681063;
    Q(3) =0.686624;
    Q(4) =0.453942;
    Q(5) =0.882871;
    return Q;
}
vector<double> get_P(){
    
    vector<double> P(6);
    
    P(0) =   -0.293855;
    P(1) =   0.179483;
    P(2) =   0.0933905;
    P(3) =   0.209977;
    P(4) =   -0.288236;
    P(5) =   -0.10935;
    return P;
}
matrix <double> get_xElec(){
    
    matrix<double> x(6,2);
    
    x(0,0) = 0.971536;
    x(0,1) =-1.7087;
    x(1,0) = -1.48713;
    x(1,1) = 1.19413;
    x(2,0) = -0.520063;
    x(2,1) = -0.715941;
    x(3,0) = 0.0667788;
    x(3,1) = -1.52122;
    x(4,0) = -0.42076;
    x(4,1) = 0.0764027;
    x(5,0) = -1.23619;
    x(5,1) = -0.022386;
    
    return x;
}
matrix<double> get_pElec(){
    
    matrix<double> p(6,2);
    
    p(0,0) =  1.0365;
    p(0,1) =  -1.32068;
    p(1,0) =  -0.797405;
    p(1,1) =  1.15937;
    p(2,0) =  0.613148;
    p(2,1) =  -0.546167;
    p(3,0) =  -1.06951;
    p(3,1) =  1.05173;
    p(4,0) =  0.240785;
    p(4,1) =  -0.7136;
    p(5,0) =  -0.640736;
    p(5,1) =  -0.163804;
    
    return p;
}
void write_v(std::ofstream &s, vector<double> v, double dt){
    
    s << dt << " ";
    
    for (int i=0; i<6; i++) {
        s << v(i) << " ";
    }
    
    s << std::endl;
}

void write_m(std::ofstream &s, matrix<double> m, double dt){
    
    s << dt << " ";
    
    for (int i=0; i<6; i++) {
        for (int j=0; j<2; j++) {
            s << m(i,j) << " ";
        }
    }
    
    s << std::endl;
}

using namespace boost::numeric::ublas;

int main(int argc, char ** argv){
    

    int nuc_beads, elec_beads, num_states;
    nuc_beads = 6;
    elec_beads = 6;
    num_states = 2;
    double beta = 1.0;
    double mass = 1.0;
    double alpha = 1.0;
    double dt = 0.1;
    double time = 10.0;

    vector<double> Q1(nuc_beads,0);
    vector<double> P1(nuc_beads,0);
    matrix<double> x1(elec_beads,num_states,0);
    matrix<double> p1(elec_beads,num_states,0);

//    vector<double> Q2(nuc_beads,0);
//    vector<double> P2(nuc_beads,0);
//    matrix<double> x2(elec_beads,num_states,0);
//    matrix<double> p2(elec_beads,num_states,0);
//
    Q1 = get_Q();
    P1 = get_P();
    x1 = get_xElec();
    p1 = get_pElec();

//    Q2 = get_Q();
//    P2 = get_P();
//    x2 = get_xElec();
//    p2 = get_pElec();
//
//    for(int bead=0; bead<nuc_beads; bead++){
//        Q1(bead) = bead - 3;
//        P1(bead) = bead*2.78;
//    }
//
//    for(int bead=0;bead<elec_beads;bead++){
//        for(int state=0;state<num_states;state++){
//            x1(bead,state) = 0.1*bead - 0.2*state;
//            p1(bead,state) = 0.1*bead + 0.3*state;
//        }
//    }
//
//
//    Q2 = Q1;
//    P2 = P1;
//    x2 = x1;
//    p2 = p1;
//
    /* New Forces*/
    C_Matrix C1(elec_beads,num_states,alpha);
    M_Matrix M1(num_states,elec_beads,beta/elec_beads);
    dM_Matrix_dQ M_dQ1(elec_beads,num_states,beta/elec_beads,M1);
    theta_mixed theta1(num_states,nuc_beads,elec_beads,C1,M1);
    theta_mixed_dQ theta_dQ1(num_states,nuc_beads,elec_beads,C1,M1,M_dQ1);
    theta_mixed_dElec theta_dElec1(num_states,elec_beads,alpha,C1,M1);
    mvrpmd_mixed_forces F1(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
                           alpha,theta1,theta_dQ1,theta_dElec1);

    ABM_MVRPMD abm1(F1,dt,num_states,nuc_beads,elec_beads);

    SpringEnergy Vspring1(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V01(nuc_beads,mass);
    GTerm G1(elec_beads,num_states,alpha);
    mvrpmd_mixed_ham H1(beta/nuc_beads,Vspring1,V01,G1,theta1);

//    /* Old Forces */
//
//    C_Matrix C2(elec_beads,num_states,alpha);
//    M_Matrix M2(num_states,nuc_beads,beta/nuc_beads);
//    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M2);
//    dM_Matrix_dQ M_dQ2(nuc_beads,num_states,beta/nuc_beads,M2);
//    dM_Matrix_MTS_dQ dM_MTS_dQ(nuc_beads,elec_beads,num_states,M_dQ2);
//    Theta_MTS theta2(num_states,elec_beads,C2,M_MTS);
//    dTheta_MTS_dQ theta_dQ2(num_states,nuc_beads,elec_beads,C2,M_MTS,dM_MTS_dQ);
//    dTheta_MTS_dElec theta_dElec2(num_states,elec_beads,C2,M_MTS);
//    Forces_MTS F2(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
//                  theta2, theta_dQ2, theta_dElec2);
//
//    ABM_MVRPMD abm2(F2,dt,num_states,nuc_beads,elec_beads);
//
//    SpringEnergy Vspring2(nuc_beads,mass,beta/nuc_beads);
//    StateIndepPot V02(nuc_beads,mass);
//    GTerm G2(elec_beads,num_states,alpha);
//    MVRPMD_MTS_Hamiltonian H2(beta/nuc_beads,Vspring2,V02,G2,theta2);
//
//
//
    int num_steps = time/dt;


    std::ofstream s_Q1, s_P1, s_p1, s_x1;
    s_Q1.open("Q1"), s_P1.open("P1"), s_x1.open("xElec1"),s_p1.open("pElec1");

    abm1.initialize_rk4(Q1,P1,x1,p1);

    F1.print_dHdQ();
    F1.print_dHdP();
    F1.print_dHdx();
    F1.print_dHdp();


    for (int i=0; i<num_steps; i++) {
        abm1.take_step(Q1,P1,x1,p1);
        write_v(s_Q1,Q1,i*dt);
        write_v(s_P1,P1,i*dt);
        write_m(s_x1,x1,i*dt);
        write_m(s_p1,p1,i*dt);
        
        if(contains_NaN(Q1)){
            std::cout << "nanned" << std::endl;
        }
    }
    //
//    s_Q1.close(), s_P1.close(), s_x1.close(),s_p1.close();
//
//
//    std::ofstream s_Q2, s_P2, s_p2, s_x2;
//    s_Q2.open("Q2"), s_P2.open("P2"), s_x2.open("xElec2"),s_p2.open("pElec2");
//
//    abm2.initialize_rk4(Q2,P2,x2,p2);
//
//    F2.print_dHdQ();
//    F2.print_dHdP();
//    F2.print_dHdx();
//    F2.print_dHdp();
//
//    for (int i=0; i<num_steps; i++) {
//        abm2.take_step(Q2,P2,x2,p2);
//        write_v(s_Q2,Q2,i*dt);
//        write_v(s_P2,P2,i*dt);
//        write_m(s_x2,x2,i*dt);
//        write_m(s_p2,p2,i*dt);
//    }
//
//    s_Q2.close(), s_P2.close(), s_x2.close(),s_p2.close();
    
    double x = -1.0;
    double y = 5.0;
    double z = x + y;
    
    std::cout << z << std::endl;
    
    std::cout << is_NaN(z) << std::endl;
    

    return 0;
}

#include "RK4_MVRPMD.hpp"

RK4_MVRPMD::RK4_MVRPMD(Forces_MTS *F_In,int nuc_beads,int elec_beads,int num_states,double dt)

    :nuc_beads(nuc_beads),num_states(num_states),elec_beads(elec_beads),
     dt(dt), dt_half(0.5*dt),
     coeff1(dt/6.0),coeff2(2.0*dt/6.0),coeff3(2.0*dt/6.0),coeff4(dt/6.0),

     k1Q(nuc_beads,0.0),k2Q(nuc_beads,0.0),k3Q(nuc_beads,0.0),k4Q(nuc_beads,0.0),
     k1P(nuc_beads,0.0),k2P(nuc_beads,0.0),k3P(nuc_beads,0.0),k4P(nuc_beads,0.0),
   
     k1x(elec_beads,num_states,0.0),k2x(elec_beads,num_states,0.0),k3x(elec_beads,num_states,0.0),
     k4x(elec_beads,num_states,0.0),

     k1p(elec_beads,num_states,0.0),k2p(elec_beads,num_states,0.0),k3p(elec_beads,num_states,0.0),
     k4p(elec_beads,num_states,0.0)

{
    F = F_In;
}

void RK4_MVRPMD::take_step(vector<double> &Q, vector<double> &P,
                           matrix<double> &x, matrix<double> &p){


    /* Round 1 */

//    std::cout << std::endl << std::endl;
//    std::cout << "RK4 1st step" << std::endl;

    F->update_Forces(Q,P,x,p);

//    F->print_dHdQ();
//    F->print_dHdP();
//    F->print_dHdp();
//    F->print_dHdx();
    
    update_k1();

    /* Round 2 */
//    std::cout << std::endl << std::endl;
//    std::cout << "RK4 2nd step" << std::endl;
//
//    std::cout << std::endl;
//    std::cout << "k1Q" << std::endl;
//    std::cout << k1Q << std::endl << std::endl;
//
//    std::cout << std::endl;
//    std::cout << "k1P" << std::endl;
//    std::cout << k1P << std::endl << std::endl;
//
//    std::cout << std::endl;
//    std::cout << "k1p" << std::endl;
//    std::cout << k1p << std::endl << std::endl;
//
//    std::cout << std::endl;
//    std::cout << "k1x" << std::endl;
//    std::cout << k1x << std::endl << std::endl;


    F->update_Forces(Q + dt_half*k1Q, P + dt_half*k1P,
                     x + dt_half*k1x, p + dt_half*k1p);
   
//    F->print_dHdQ();
//    F->print_dHdP();
//    F->print_dHdp();
//    F->print_dHdx();

    update_k2();
  
  
    /* Round 3 */
//    std::cout << std::endl << std::endl;
//    std::cout << "RK4 3rd step" << std::endl;
    
    F->update_Forces(Q + dt_half*k2Q, P + dt_half*k2P,
                     x + dt_half*k2x, p + dt_half*k2p);
 
//    F->print_dHdQ();
//    F->print_dHdP();
//    F->print_dHdp();
//    F->print_dHdx();
    
    update_k3();
    

    /* Round 4 */

//    std::cout << std::endl << std::endl;
//    std::cout << "RK4 4th step" << std::endl;
    
//    F->print_dHdQ();
//    F->print_dHdP();
//    F->print_dHdp();
//    F->print_dHdx();
    
    F->update_Forces(Q + dt*k3Q, P + dt*k3P,
                     x + dt*k3x, p + dt*k3p );
    update_k4();
    
    update_final(Q,P,x,p);
}

void RK4_MVRPMD::update_k1(){
    
    k1Q = F->get_dHdP();
    k1P = - F->get_dHdQ();
    k1x =  F->get_dHdp();
    k1p = - F->get_dHdx();

}

void RK4_MVRPMD::update_k2(){
    
    k2Q = F->get_dHdP();
    k2P = - F->get_dHdQ();
    k2x =  F->get_dHdp();
    k2p = - F->get_dHdx();
}

void RK4_MVRPMD::update_k3(){
    
    k3Q = F->get_dHdP();
    k3P = - F->get_dHdQ();
    k3x =  F->get_dHdp();
    k3p = - F->get_dHdx();
}

void RK4_MVRPMD::update_k4(){
    
    k4Q = F->get_dHdP();
    k4P = - F->get_dHdQ();
    k4x =  F->get_dHdp();
    k4p = - F->get_dHdx();
}

void RK4_MVRPMD::update_final(vector<double> &Q, vector<double> &P,
                              matrix<double> &x,matrix<double> &p){

  Q = Q + (coeff1 * k1Q) + (coeff2 * k2Q) + (coeff3 * k3Q) + (coeff4 *k4Q);
  P = P + (coeff1 * k1P) + (coeff2 * k2P) + (coeff3 * k3P) + (coeff4 *k4P);
  x  = x  + (coeff1 * k1x)  + (coeff2 * k2x)  + (coeff3 * k3x)  + (coeff4 *k4x);
  p  = p  + (coeff1 * k1p)  + (coeff2 * k2p)  + (coeff3 * k3p)  + (coeff4 *k4p);
}

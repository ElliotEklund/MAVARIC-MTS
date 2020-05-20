#include "ABM_MVRPMD.hpp"
#include <iomanip>


ABM_MVRPMD::ABM_MVRPMD(Forces_MTS &F_In, double dt, int num_states, int nuc_beads, int elec_beads)
    :dt(dt), num_states(num_states), elec_beads(elec_beads), nuc_beads(nuc_beads),
     h1_p(dt*55.0/24.0),h2_p(dt*59.0/24.0),h3_p(dt*37.0/24.0),h4_p(dt*9.0/24.0),
     h1_c(dt*9.0/24.0),h2_c(dt*19.0/24.0),h3_c(dt*5.0/24.0),h4_c(dt*1.0/24.0),

     f_Q_3(nuc_beads,0.0),f_Q_2(nuc_beads,0.0),f_Q_1(nuc_beads,0.0),f_Q_0(nuc_beads,0.0),
     f_Q_p1(nuc_beads,0.0),

     f_P_3(nuc_beads,0.0),f_P_2(nuc_beads,0.0),f_P_1(nuc_beads,0.0),f_P_0(nuc_beads,0.0),
     f_P_p1(nuc_beads,0.0),

     f_x_3(elec_beads,num_states,0.0),f_x_2(elec_beads,num_states,0.0),f_x_1(elec_beads,num_states,0.0),
     f_x_0(elec_beads,num_states,0.0),f_x_p1(elec_beads,num_states,0.0),

     f_p_3(elec_beads,num_states,0.0),f_p_2(elec_beads,num_states,0.0),f_p_1(elec_beads,num_states,0.0),
     f_p_0(elec_beads,num_states,0.0),f_p_p1(elec_beads,num_states,0.0),

     pred_Q(nuc_beads,0.0),pred_P(nuc_beads,0.0),pred_x(elec_beads,num_states,0.0),pred_p(elec_beads,num_states,0.0)

{
    F = &F_In;
}

 void ABM_MVRPMD::initialize_rk4(vector<double> &Q,vector<double> &P,
   matrix<double> &x,matrix<double> &p){

  vector<double> Q_copy = Q;
  vector<double> P_copy = P;
  matrix<double> x_copy = x;
  matrix<double> p_copy = p;

  /* Use RK4 to take steps backwards. These steps will seed ABM. */
  RK4_MVRPMD my_rk4(F,nuc_beads,elec_beads,num_states,-dt);

//  std::cout << std::endl << std::endl;
//  std::cout << "-----------------------" << std::endl;
//  std::cout << "1st step back" << std::endl;

  my_rk4.take_step(Q,P,x,p);
  update_f_1(Q,P,x,p);

  /* Take step backward to t_(-1)*/

//  std::cout << std::endl << std::endl;
//  std::cout << "-----------------------" << std::endl;
//  std::cout << "2nd step back" << std::endl;
  
  my_rk4.take_step(Q,P,x,p);
  update_f_2(Q,P,x,p);

  /* Take step backward to t_(-2)*/

//  std::cout << std::endl << std::endl;
//  std::cout << "-----------------------" << std::endl;
//  std::cout << "3rd step back" << std::endl;
  
  my_rk4.take_step(Q,P,x,p);
  update_f_3(Q,P,x,p);

  /* Copy original PSV at t_0 back to outgoing variables. */
  noalias(Q) = Q_copy;
  noalias(P) = P_copy;
  noalias(x) = x_copy;
  noalias(p) = p_copy;
}


void ABM_MVRPMD::take_step(vector<double> &Q, vector<double> &P,
                           matrix<double> &x, matrix<double> &p){
    
//    std::cout << "Predict step " << std::endl << std::endl;
    
    predict(Q,P,x,p);
 
//    std::cout << "Q predict" << std::endl;
//    std::cout << pred_Q << std::endl << std::endl;
//
//    std::cout << "P predict" << std::endl;
//    std::cout << pred_P << std::endl << std::endl;
//
//    std::cout << "x predict" << std::endl;
//    std::cout << pred_x << std::endl << std::endl;
//
//    std::cout << "p predict" << std::endl;
//    std::cout << pred_p << std::endl << std::endl;
//
//    std::cout << "Correct step " << std::endl << std::endl;
    
    correct(Q,P,x,p);

}


void ABM_MVRPMD::predict(const vector<double> &Q,const vector<double> &P,
                         const matrix<double> &x,const matrix<double> &p){
    
    update_f_0(Q,P,x,p);
    
    noalias(pred_Q) = Q + (h1_p * f_Q_0) - (h2_p * f_Q_1) + (h3_p * f_Q_2) - (h4_p * f_Q_3);
    noalias(pred_P) = P + (h1_p * f_P_0) - (h2_p * f_P_1) + (h3_p * f_P_2) - (h4_p * f_P_3);
    noalias(pred_x) = x + (h1_p * f_x_0) - (h2_p * f_x_1) + (h3_p * f_x_2) - (h4_p * f_x_3);
    noalias(pred_p) = p + (h1_p * f_p_0) - (h2_p * f_p_1) + (h3_p * f_p_2) - (h4_p * f_p_3);
    
//    pred_Q = Q + (h1_p * f_Q_0) - (h2_p * f_Q_1) + (h3_p * f_Q_2) - (h4_p * f_Q_3);
//    pred_P = P + (h1_p * f_P_0) - (h2_p * f_P_1) + (h3_p * f_P_2) - (h4_p * f_P_3);
//    pred_x = x + (h1_p * f_x_0) - (h2_p * f_x_1) + (h3_p * f_x_2) - (h4_p * f_x_3);
//    pred_p = p + (h1_p * f_p_0) - (h2_p * f_p_1) + (h3_p * f_p_2) - (h4_p * f_p_3);
    
//    std::cout << "p in predict" << std::endl;
//    std::cout << p << std::endl;
//    
//    std::cout << "f_p_0" << std::endl;
//    std::cout << h1_p * f_p_0 << std::endl;
//    
//    std::cout << "f_p_1" << std::endl;
//    std::cout << h2_p * f_p_1 << std::endl;
//    
//    std::cout << "f_p_2" << std::endl;
//    std::cout << h3_p * f_p_2 << std::endl;
    
//    std::cout << std::scientific << std::setprecision(15) << "f_p_3" << std::endl;
//    std::cout << h4_p * f_p_3 << std::endl;

}

void ABM_MVRPMD::correct(vector<double> &Q, vector<double> &P,
                         matrix<double> &x, matrix<double> &p){
    
    update_f_p1(pred_Q,pred_P,pred_x,pred_p);
    
    Q  = Q  + (h1_c * f_Q_p1)  + (h2_c * f_Q_0)  - (h3_c * f_Q_1)  + (h4_c * f_Q_2);
    P  = P  + (h1_c * f_P_p1)  + (h2_c * f_P_0)  - (h3_c * f_P_1)  + (h4_c * f_P_2);
    x  = x  + (h1_c * f_x_p1)  + (h2_c * f_x_0)  - (h3_c * f_x_1)  + (h4_c * f_x_2);
    p  = p  + (h1_c * f_p_p1)  + (h2_c * f_p_0)  - (h3_c * f_p_1)  + (h4_c * f_p_2);
    
    noalias(f_Q_3) = f_Q_2;
    noalias(f_Q_2) = f_Q_1;
    noalias(f_Q_1) = f_Q_0;

    noalias(f_P_3) = f_P_2;
    noalias(f_P_2) = f_P_1;
    noalias(f_P_1) = f_P_0;
  
    noalias(f_x_3) = f_x_2;
    noalias(f_x_2) = f_x_1;
    noalias(f_x_1) = f_x_0;

    noalias(f_p_3) = f_p_2;
    noalias(f_p_2) = f_p_1;
    noalias(f_p_1) = f_p_0;
    
}

void ABM_MVRPMD::update_f_3(const vector<double> &Q,const vector<double> &P,const matrix<double> &x,
                            const matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    
    noalias(f_Q_3)  =   F->get_dHdP();
    noalias(f_P_3)  = - F->get_dHdQ();
    noalias(f_x_3)  =   F->get_dHdp();
    noalias(f_p_3)  = - F->get_dHdx();
}

void ABM_MVRPMD::update_f_2(const vector<double> &Q,const vector<double> &P,const matrix<double> &x,
                            const matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    
    noalias(f_Q_2)  =   F->get_dHdP();
    noalias(f_P_2)  = - F->get_dHdQ();
    noalias(f_x_2)  =   F->get_dHdp();
    noalias(f_p_2)  = - F->get_dHdx();
}

void ABM_MVRPMD::update_f_1(const vector<double> &Q,const vector<double> &P,const matrix<double> &x,
                            const matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    
    noalias(f_Q_1) =   F->get_dHdP();
    noalias(f_P_1) = - F->get_dHdQ();
    noalias(f_x_1)  =   F->get_dHdp();
    noalias(f_p_1)  = - F->get_dHdx();
}

void ABM_MVRPMD::update_f_0(const vector<double> &Q,const vector<double> &P,const matrix<double> &x,
                            const matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    
    noalias(f_Q_0) =  F->get_dHdP();
    noalias(f_P_0) = -F->get_dHdQ();
    noalias(f_x_0) =  F->get_dHdp();
    noalias(f_p_0) = -F->get_dHdx();
}

void ABM_MVRPMD::update_f_p1(const vector<double> &Q,const vector<double> &P,const matrix<double> &x,
                             const matrix<double> &p){
    
    F->update_Forces(Q,P,x,p);
    
    noalias(f_Q_p1) =  F->get_dHdP();
    noalias(f_P_p1) = -F->get_dHdQ();
    noalias(f_x_p1) =  F->get_dHdp();
    noalias(f_p_p1) = -F->get_dHdx();
}

#include "Theta.h"

Theta::Theta(int num_states, int num_beads, double beta_num_beads,
             C_Matrix &C_In, M_Matrix &M_In)
    :num_states(num_states), num_beads(num_beads), beta_num_beads(beta_num_beads),
     gamma_mat(num_states,num_states),
     temp1(num_states,num_states),temp2(num_states,num_states),myI(num_states,num_states)
{
    C = &C_In;
    M = &M_In;
    gamma_mat = identity_matrix<std::complex<double> > (num_states);
    myI = identity_matrix<std::complex<double> > (num_states);
}

void Theta::update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                             const matrix<double> &p){
    
    C->update_C_vec(x, p);
    M->update_M_vec(Q);

    temp1 = myI;
    temp2 = myI;

    /* Chain multiplicaton */
    for(int bead=0; bead<num_beads; bead++){
        
       // gamma_mat = prod(gamma_mat,C->get_C_alpha(bead));
       // gamma_mat = prod(gamma_mat,M->get_M_alpha(bead));
        noalias(temp2) = prod(temp1,C->get_C_alpha(bead));
        noalias(temp1) = prod(temp2,M->get_M_alpha(bead));
    }

    noalias(gamma_mat) = temp1;
}

void Theta::update_theta(const vector<double> &Q,const matrix<double> &x,
                         const matrix<double> &p){
    
    update_gamma_mat(Q, x, p);
    
    /* Compute trace. */
    std::complex<double> tr(0,0);
    for(int state=0; state<num_states; state++){
      tr += gamma_mat(state,state);
    }
    
    theta = tr.real();
}

double& Theta::get_theta(const vector<double> &Q,const matrix<double> &x,
                        const matrix<double> &p){
    
    update_theta(Q, x, p);
    return theta;
}

double& Theta::get_theta(){return theta;}

double Theta::get_signTheta(){
    if (theta >= 0) {
        return 1.0;
    }
    else{
        return -1.0;
    }
}

#include "dTheta_dBeta.hpp"

dTheta_dBeta::dTheta_dBeta(int num_beads, int num_states, double beta_num_beads,
                           C_Matrix &C_In, M_Matrix &M_In)
    :dM_dBeta(num_beads,num_states,beta_num_beads,M_In),
     f_chain(num_beads), b_chain(num_beads), num_beads(num_beads),
     num_states(num_states), dGamma_dBeta(num_states,num_states),
     f_temp(num_states,num_states),b_temp(num_states,num_states),
     dummy1(num_states,num_states),dummy2(num_states,num_states)
{
    C = &C_In;
    M = &M_In;
}

double dTheta_dBeta::get_dTheta_dBeta(){
    update_dTheta_dBeta();
    return d_Theta_dBeta;
}

void dTheta_dBeta::update_dTheta_dBeta(){

    update_dGamma_dBeta();
    std::complex<double> tr (0,0); //compute trace

    for (int state=0; state<num_states; state++) {
        tr += dGamma_dBeta(state,state);
    }

    /* Compute trace. */
    d_Theta_dBeta = tr.real();
}

void dTheta_dBeta::update_dGamma_dBeta(){
    
    update_f_chain();
    update_b_chain();
    dM_dBeta.update_dM_dBeta_vec();

    matrix<std::complex<double> > sum (num_states,num_states);
    std::fill(sum.data().begin(),sum.data().end(),0.0);
    
    for (int bead=0; bead<num_beads; bead++) {
        noalias(dummy1) = prod(f_chain(bead), dM_dBeta.get_dM_dBeta_alpha(bead));
        noalias(dummy2) = prod(dummy1, b_chain(bead));
        sum += dummy2;
    }
    noalias(dGamma_dBeta) = sum;
}

void dTheta_dBeta::update_f_chain(){

    f_chain[0] = C->get_C_alpha(0);
    
    for (int bead=1; bead<num_beads; bead++) {
        noalias(f_temp) = prod(M->get_M_alpha(bead-1), C->get_C_alpha(bead));
        f_chain[bead] = prod(f_chain[bead-1], f_temp);
    }
}

void dTheta_dBeta::update_b_chain(){
    
    b_chain[num_beads-1] = identity_matrix<std::complex<double> > (num_states);
    
    for (int bead = num_beads - 2; bead > -1; bead--) {
        noalias(b_temp) = prod(C->get_C_alpha(bead+1), M->get_M_alpha(bead+1));
        b_chain[bead] = prod(b_temp, b_chain[bead+1]);
    }
}

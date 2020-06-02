#include "dTheta_MTS_dBeta.hpp"

dTheta_MTS_dBeta::dTheta_MTS_dBeta(int nuc_beads, int elec_beads, int num_states,
                                   double beta_num_beads,C_Matrix &C_In,
                                   M_Matrix &M_In, M_Matrix_MTS &M_MTS_In)
    :dM_MTS_dBeta(nuc_beads,elec_beads,num_states,beta_num_beads,M_In),
     f_chain(elec_beads), b_chain(elec_beads), nuc_beads(nuc_beads),
     elec_beads(elec_beads),num_states(num_states),
     dGamma_dBeta(num_states,num_states)
{
    C = &C_In;
    M_MTS = &M_MTS_In;
}
double dTheta_MTS_dBeta::get_dTheta_MTS_dBeta(){
    update_dTheta_MTS_dBeta();
    return d_Theta_MTS_dBeta;
}
void dTheta_MTS_dBeta::update_dTheta_MTS_dBeta(){

    update_dGamma_dBeta();
    std::complex<double> tr (0,0); //compute trace

    for (int state=0; state<num_states; state++) {
        tr += dGamma_dBeta(state,state);
    }

    /* Compute trace. */
    d_Theta_MTS_dBeta = tr.real();
}
void dTheta_MTS_dBeta::update_dGamma_dBeta(){
    
    update_f_chain();
    update_b_chain();
    
    dM_MTS_dBeta.update_dM_MTS_dBeta_vec();

    matrix<std::complex<double> > temp (num_states,num_states);
    matrix<std::complex<double> > sum (num_states,num_states);
    
    std::fill(sum.data().begin(),sum.data().end(),0.0);
    
    for (int bead=0; bead<elec_beads; bead++) {
        temp = prod(f_chain(bead), dM_MTS_dBeta.get_dM_MTS_dBeta_alpha(bead));
        temp = prod(temp, b_chain(bead));
        sum += temp;
    }
    dGamma_dBeta = sum;
}
void dTheta_MTS_dBeta::update_f_chain(){
    
    f_chain[0] = C->get_C_alpha(0);
    matrix<std::complex<double> > temp(num_states,num_states);
    
    for (int bead=1; bead<elec_beads; bead++) {
        temp = prod(M_MTS->get_M_MTS_alpha(bead-1), C->get_C_alpha(bead));
        f_chain[bead] = prod(f_chain[bead-1], temp);
    }
}
void dTheta_MTS_dBeta::update_b_chain(){
    
    b_chain[elec_beads-1] = identity_matrix<std::complex<double> > (num_states);
    matrix<std::complex<double> > temp(num_states,num_states);
    
    for (int bead = elec_beads - 2; bead > -1; bead--) {
        temp = prod(C->get_C_alpha(bead+1), M_MTS->get_M_MTS_alpha(bead+1));
        b_chain[bead] = prod(temp, b_chain[bead+1]);
    }
}

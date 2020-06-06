#include "theta_mixed_dBeta.hpp"

theta_mixed_dBeta::theta_mixed_dBeta(int elec_beads, int num_states,
                                   double beta_num_beads,C_Matrix &C_In,
                                   M_Matrix &M_In)
    :f_chain(elec_beads), b_chain(elec_beads),
     elec_beads(elec_beads),num_states(num_states),
     dgamma(num_states,num_states),
     M_dbeta(elec_beads,num_states,beta_num_beads,M_In)
{
    C = &C_In;
    M = &M_In;
}
double theta_mixed_dBeta::get_dtheta(){
    update_dtheta();
    return dtheta;
}
void theta_mixed_dBeta::update_dtheta(){

    update_dgamma();
    std::complex<double> tr (0,0);
    tr = trace<std::complex<double> >(dgamma,num_states);
    dtheta = tr.real();
}
void theta_mixed_dBeta::update_dgamma(){

    update_f_chain();
    update_b_chain();

    M_dbeta.update_dM_dBeta_vec();

    matrix<std::complex<double> > temp (num_states,num_states,0);
    matrix<std::complex<double> > sum (num_states,num_states,0);

    for (int bead=0; bead<elec_beads; bead++) {
        temp = prod(f_chain(bead), M_dbeta.get_dM_dBeta_alpha(bead));
        temp = prod(temp, b_chain(bead));
        sum += temp;
    }
    dgamma = sum;
}
void theta_mixed_dBeta::update_f_chain(){

    f_chain[0] = C->get_C_alpha(0);
    matrix<std::complex<double> > temp(num_states,num_states);

    for (int bead=1; bead<elec_beads; bead++) {
        temp = prod(M->get_M_alpha(bead-1), C->get_C_alpha(bead));
        f_chain[bead] = prod(f_chain[bead-1], temp);
    }
}
void theta_mixed_dBeta::update_b_chain(){

    b_chain[elec_beads-1] = identity_matrix<std::complex<double> > (num_states);
    matrix<std::complex<double> > temp(num_states,num_states);

    for (int bead = elec_beads - 2; bead > -1; bead--) {
        temp = prod(C->get_C_alpha(bead+1), M->get_M_alpha(bead+1));
        b_chain[bead] = prod(temp, b_chain[bead+1]);
    }
}

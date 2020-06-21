#include "theta_mixed_dElec.hpp"

theta_mixed_dElec::theta_mixed_dElec(int num_states, int elec_beads, double alpha,
                                     C_Matrix &C_In, M_Matrix &M_In)
    :num_states(num_states),elec_beads(elec_beads),
     dC(elec_beads, num_states,alpha),

     f_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),
     b_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),

     theta_dx_vec(elec_beads,num_states,0.0),
     theta_dp_vec(elec_beads,num_states,0.0),

     xtemp1(num_states,num_states,0.0),xtemp2(num_states,num_states,0.0),
     ptemp1(num_states,num_states,0.0),ptemp2(num_states,num_states,0.0),

     fchain_temp(num_states,num_states,0.0),bchain_temp(num_states,num_states,0.0)
{
    C = &C_In;
    M = &M_In;
}
void theta_mixed_dElec::update_theta_dElec(const matrix<std::complex<double> > & x,
                                           const matrix<std::complex<double> > & p){
    update_f_chain();
    update_b_chain();
    update_theta_dx_vec(x,p);
    update_theta_dp_vec(x,p);
}
void theta_mixed_dElec::update_theta_dx_vec(const matrix<std::complex<double> > & x,
                                            const matrix<std::complex<double> > & p){
    dC.update_dCdx_mat(x,p);
    std::complex<double> tr (0,0); //compute trace
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            noalias(xtemp1) = prod(f_chain(bead), dC.get_dC_dx(bead,state));
            noalias(xtemp2) = prod(xtemp1, b_chain(bead));
            tr = trace<std::complex<double> >(xtemp2,num_states);
            theta_dx_vec(bead,state) = tr.real();
            tr = std::complex<double> (0,0);
        }
    }
}
void theta_mixed_dElec::update_theta_dp_vec(const matrix<std::complex<double> > & x,
                                            const matrix<std::complex<double> > & p){
    
    dC.update_dCdp_mat(x,p);
    std::complex<double> tr (0,0); //compute trace
        
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            noalias(ptemp1) = prod(f_chain(bead), dC.get_dC_dp(bead,state));
            noalias(ptemp2) = prod(ptemp1, b_chain(bead));
            tr = trace<std::complex<double> >(ptemp2,num_states);
            theta_dp_vec(bead,state) = tr.real();
            tr = std::complex<double> (0,0);
        }
    }
}
void theta_mixed_dElec::update_f_chain(){
    
    f_chain[0] = identity_matrix<std::complex<double> > (num_states);
    
    for (int bead=1; bead<elec_beads; bead++) {
        noalias(fchain_temp) = prod(C->get_C_alpha(bead-1),M->get_M_alpha(bead-1));
        f_chain[bead] = prod(f_chain[bead-1], fchain_temp);
    }
}
void theta_mixed_dElec::update_b_chain(){
    
    b_chain[elec_beads-1] = M->get_M_alpha(elec_beads-1);
    
    for (int bead = elec_beads - 2; bead > -1; bead--) {
        noalias(bchain_temp) = prod(M->get_M_alpha(bead),C->get_C_alpha(bead+1));
        b_chain[bead] = prod(bchain_temp, b_chain[bead+1]);
    }
}
const matrix<double> & theta_mixed_dElec::get_theta_dx_vec(){
    return theta_dx_vec;}
const matrix<double> & theta_mixed_dElec::get_theta_dp_vec(){
    return theta_dp_vec;}

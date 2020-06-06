#include "dTheta_MTS_dQ.hpp"

dTheta_MTS_dQ::dTheta_MTS_dQ(int num_states, int nuc_beads, int elec_beads, C_Matrix &C_In,
                             M_Matrix_MTS &M_MTS_In,dM_Matrix_MTS_dQ &dM_MTS_dQ_In)

    :num_states(num_states), elec_beads(elec_beads),
     nuc_beads(nuc_beads), gamma(nuc_beads/elec_beads),

     f_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),
     b_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),

     dTheta_MTS_dQ_vec(nuc_beads,0.0),

     temp1(num_states,num_states,0.0),temp2(num_states,num_states,0.0),
     fchain_temp(num_states,num_states,0.0),bchain_temp(num_states,num_states,0.0)

{
    C = &C_In;
    M_MTS = &M_MTS_In;
    dM_MTS_dQ = & dM_MTS_dQ_In;
}

void dTheta_MTS_dQ::update_dTheta_MTS_dQ_vec(const vector<double> &Q){
      
    update_f_chain();
    update_b_chain();

    dM_MTS_dQ->update_dM_MTS_dQ_vec(Q);
    int alpha = 0;
    
    std::complex<double> tr (0,0); //compute trace
    /* It is assummed that dM_MTS_dQ has already been updated*/
    
    for (int bead=0; bead<nuc_beads; bead++) {
        
        alpha = floor(bead/gamma);
        
        noalias(temp1) = prod(f_chain(alpha), dM_MTS_dQ->get_dM_MTS_dQ_alpha(bead));
        noalias(temp2) = prod(temp1, b_chain(alpha));

        for (int state=0; state<num_states; state++) {
            tr += temp2(state,state);
        }
        
        /* Compute trace. */
        dTheta_MTS_dQ_vec(bead) = tr.real();
        tr = std::complex<double> (0,0);
    }
}
void dTheta_MTS_dQ::update_f_chain(){
    
    f_chain[0] = C->get_C_alpha(0);
    
    for (int bead=1; bead<elec_beads; bead++) {
        noalias(fchain_temp) = prod(M_MTS->get_M_MTS_alpha(bead-1), C->get_C_alpha(bead));
        f_chain[bead] = prod(f_chain[bead-1], fchain_temp);
    }
}
void dTheta_MTS_dQ::update_b_chain(){
    
    b_chain[elec_beads-1] = identity_matrix<std::complex<double> > (num_states);
    
    for (int bead = elec_beads - 2; bead > -1; bead--) {
        noalias(bchain_temp) = prod(C->get_C_alpha(bead+1), M_MTS->get_M_MTS_alpha(bead+1));
        b_chain[bead] = prod(bchain_temp, b_chain[bead+1]);
    }
}
const vector<double> & dTheta_MTS_dQ::get_dThetaMTS_dQ_vec(){
    return dTheta_MTS_dQ_vec;
}

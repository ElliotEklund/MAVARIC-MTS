#include "dM_Matrix_MTS_dQ.hpp"

dM_Matrix_MTS_dQ::dM_Matrix_MTS_dQ(int nuc_beads, int elec_beads, int num_states,
                                   dM_Matrix_dQ &dMdQ_In)
    :nuc_beads(nuc_beads),elec_beads(elec_beads),
     gamma(nuc_beads/elec_beads), ONE_gamma(1.0/gamma),
     dM_MTS_dQ_vec(nuc_beads,zero_matrix<std::complex<double> >(num_states,num_states))
{
    dMdQ = &dMdQ_In;
}
void dM_Matrix_MTS_dQ::update_dM_MTS_dQ_vec(const vector<double> &Q){
    dMdQ->update_dM_dQ_vec(Q);
}
matrix<std::complex<double> > dM_Matrix_MTS_dQ::get_dM_MTS_dQ_alpha(int alpha){
    return ONE_gamma * dMdQ->get_dMdQ_alpha(alpha);
}

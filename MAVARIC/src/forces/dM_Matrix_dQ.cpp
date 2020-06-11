#include "dM_Matrix_dQ.hpp"

dM_Matrix_dQ::dM_Matrix_dQ(int nuc_beads,int num_states, double betaN,
                           M_Matrix &M_In)
    :nuc_beads(nuc_beads), num_states(num_states), betaN(betaN),

     dVdQ(nuc_beads, num_states),

     dVdQ_mat(nuc_beads,num_states,0.0),
     dM_dQ_alpha(num_states,num_states,0.0),
     dVcouple_dQ_vec(nuc_beads,zero_matrix<double>(num_states,num_states)),
     dM_dQ_vec(nuc_beads,zero_matrix<std::complex<double> >(num_states,num_states))
{
    M = &M_In;
}
void dM_Matrix_dQ::update_dM_dQ_vec(const vector<double> &Q){
    
    dVdQ.update_dVdQ(Q);
    dVdQ_mat = dVdQ.get_dVdQ();
    dVcouple_dQ_vec = dVdQ.get_dVcoup_dQ_vec(Q);
        
    for (int bead=0; bead<nuc_beads; bead++) {
        update_dM_dQ_alpha(Q(bead), bead);
        dM_dQ_vec(bead) = dM_dQ_alpha;
    }
}
void dM_Matrix_dQ::update_dM_dQ_alpha(const double &Q, int alpha){
    
    std::complex<double> m_diag (0,0);
    std::complex<double> term1 (0,0);
    std::complex<double> term2 (0,0);
            
    /* Update diagonal */
    for (int state=0; state<num_states; state++) {
        m_diag =  M->get_M_alpha(alpha)(state,state);
        dM_dQ_alpha(state,state) = -betaN * dVdQ_mat(alpha,state) * m_diag;
    }
    
   /* Update Top Right Triangle */
    for (int state1 = 0; state1 <num_states - 1; state1++) {
        for (int state2 = 1 + state1; state2<num_states; state2++) {
            term1 = dVcouple_dQ_vec(alpha)(state1,state2) * M->get_M_alpha(alpha)(state1,state2);
            term2 = dM_dQ_alpha(state1,state1) * M->get_V_couple_alpha(alpha)(state1,state2);
            
            dM_dQ_alpha(state1,state2) = -betaN*(term1 + term2);
        }
    }
    
    /* Update Lower Left Triangle*/
    for (int state1=1; state1<num_states; state1++) {
        for (int state2=0; state2<state1; state2++) {
            term1 = dVcouple_dQ_vec(alpha)(state1,state2) * M->get_M_alpha(alpha)(state1,state2);
            term2 = dM_dQ_alpha(state1,state1) * M->get_V_couple_alpha(alpha)(state1,state2);
            
            dM_dQ_alpha(state1,state2) = -betaN*(term1 + term2);
        }
    }
}
const vector<matrix<std::complex<double> > > & dM_Matrix_dQ::get_dMdQ_vec(){
    return dM_dQ_vec;
}
const matrix<std::complex<double> > & dM_Matrix_dQ::get_dMdQ_alpha(int alpha){
    return dM_dQ_vec(alpha);
}

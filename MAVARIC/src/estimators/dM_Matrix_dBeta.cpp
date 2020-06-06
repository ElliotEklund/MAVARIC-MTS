#include "dM_Matrix_dBeta.hpp"

dM_Matrix_dBeta::dM_Matrix_dBeta(int nuc_beads, int num_states, double betaN,
                                 M_Matrix &M_In)
    :num_states(num_states), nuc_beads(nuc_beads),
     betaN(betaN), ONE_N(1.0/nuc_beads),ONE_beta(1.0/(betaN*nuc_beads)),
     dM_dBeta_alpha(num_states,num_states),
     dM_dBeta_vec(nuc_beads)
{
    M = &M_In;
    dM_dBeta_alpha = identity_matrix<std::complex<double> > (num_states,num_states);
}
dM_Matrix_dBeta::dM_Matrix_dBeta(int nuc_beads, int elec_beads, int num_states, double
                                 betaN,M_Matrix &M_In)
    :num_states(num_states), nuc_beads(nuc_beads), elec_beads(elec_beads),
     betaN(betaN), ONE_N(1.0/elec_beads),ONE_beta(1.0/(betaN*elec_beads)),
     dM_dBeta_alpha(num_states,num_states),
     dM_dBeta_vec(nuc_beads)
{
    M = &M_In;
    dM_dBeta_alpha = identity_matrix<std::complex<double> > (num_states,num_states);
}
void dM_Matrix_dBeta::update_dM_dBeta(int alpha){
    
    std::complex<double> diag (0,0);
    std::complex<double> term1 (0,0);
    std::complex<double> term2 (0,0);
    double Vij;

    /* Update diagonal elements*/
    for (int state=0; state<num_states; state++) {
        diag = M->get_M_alpha(alpha)(state,state);
        diag = - ONE_N * diag * M->get_V_mat()(alpha,state);
        dM_dBeta_alpha(state,state) = diag;
    }

    /* Update upper right triangle */
    for (int state1 = 0; state1 <num_states - 1; state1++) {
        term1 = M->get_M_alpha(alpha)(state1,state1); //divide by beta
        term2 = dM_dBeta_alpha(state1,state1);

        for (int state2 = 1 + state1; state2<num_states; state2++) {
            Vij = M->get_V_couple_alpha(alpha)(state1,state2);
            dM_dBeta_alpha(state1,state2) = -Vij*(ONE_N*term1 + betaN*term2);
        }
    }
    
    /* Update lower left triangle*/
    for (int state1=1; state1<num_states; state1++) {
        term1 = M->get_M_alpha(alpha)(state1,state1); //divide by beta
        term2 = dM_dBeta_alpha(state1,state1); //gt vcouple
        
        for (int state2=0; state2<state1; state2++) {
            Vij = M->get_V_couple_alpha(alpha)(state1,state2);
            dM_dBeta_alpha(state1,state2) = -Vij*(ONE_N*term1 + betaN*term2);
        }
    }
}
const matrix<std::complex<double> > & dM_Matrix_dBeta::get_dM_dBeta_alpha (int alpha){
    return dM_dBeta_vec(alpha);
}
void dM_Matrix_dBeta::update_dM_dBeta_vec(){
    for (int bead=0; bead<nuc_beads; bead++) {
        update_dM_dBeta(bead);
        dM_dBeta_vec(bead) = dM_dBeta_alpha;
    }
}

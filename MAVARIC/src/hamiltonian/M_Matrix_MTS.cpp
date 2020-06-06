#include "M_Matrix_MTS.hpp"

M_Matrix_MTS::M_Matrix_MTS(int nuc_beads,int elec_beads,int num_states,
                           M_Matrix & M_In)
    :nuc_beads(nuc_beads),elec_beads(elec_beads),num_states(num_states),
     gamma(nuc_beads/elec_beads),

     M_MTS_vec(elec_beads,matrix<std::complex<double> > (num_states,num_states,0.0))
{
    M = &M_In;
}
void M_Matrix_MTS::update_M_MTS_vec(const vector<double> &Q){
    
    M->update_M_vec(Q);
    matrix<std::complex<double> > temp_sum (num_states,num_states,0.0);
    
    for (int bead=0; bead<elec_beads; bead++) {
        std::fill(temp_sum.data().begin(), temp_sum.data().end(), 0);
        for (int slice=0; slice<gamma; slice++) {
            temp_sum += M->get_M_alpha(bead*gamma + slice);
        }
        M_MTS_vec(bead) = (1.0/gamma) * temp_sum;
    }
}
void M_Matrix_MTS::update_M_MTS_vec(){
    
    matrix<std::complex<double> > temp_sum (num_states,num_states,0.0);
    
    for (int bead=0; bead<elec_beads; bead++) {
        std::fill(temp_sum.data().begin(), temp_sum.data().end(), 0);
        for (int slice=0; slice<gamma; slice++) {
            temp_sum += M->get_M_alpha(bead*gamma + slice);
        }
        M_MTS_vec(bead) = (1.0/gamma) * temp_sum;
    }
}

matrix<std::complex<double> >& M_Matrix_MTS::get_M_MTS_alpha(int alpha){
    return M_MTS_vec(alpha);
}

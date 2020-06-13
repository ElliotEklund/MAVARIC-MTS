#include "C_Matrix.h"

C_Matrix::C_Matrix(int num_beads, int num_states, double alpha)
    :num_beads(num_beads), num_states(num_states),alpha(alpha),
     unit_complex(0.0,1.0), half_identity(num_states,num_states),

    x_plus_ip_mat(num_beads,num_states,0.0),
    x_min_ip_mat(num_beads,num_states,0.0),
    x_plus_ip_mat_trans(num_states,num_beads,0.0),

    x_mat_c(num_beads,num_states,0.0),
    p_mat_c(num_beads,num_states,0.0),

    C_vec(num_beads,identity_matrix<std::complex<double> >(num_states))
{
    half_identity = 0.5 * identity_matrix<std::complex<double> > (num_states);
    x_alpha = sqrt(alpha);
    p_alpha = 1.0/sqrt(alpha);
}
void C_Matrix::update_C_vec(const matrix<double> &x_mat, const matrix<double> &p_mat){
    
    //copy x-mat,p-mat to complex matricies
    x_mat_c = x_alpha*x_mat;
    p_mat_c = p_alpha*p_mat;

    noalias(x_plus_ip_mat) = x_mat_c + unit_complex*p_mat_c;
    noalias(x_min_ip_mat) = x_mat_c - unit_complex*p_mat_c;

    noalias(x_plus_ip_mat_trans) = trans(x_plus_ip_mat);

    for (int bead=0; bead<num_beads; bead++) {
        noalias(C_vec(bead)) = outer_prod(column(x_plus_ip_mat_trans, bead), row(x_min_ip_mat, bead)) - half_identity;
    }
}
matrix<std::complex<double> >& C_Matrix::get_C_alpha(int alpha){
    return C_vec(alpha);
}

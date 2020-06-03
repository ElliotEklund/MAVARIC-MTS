#include "M_Matrix.h"

M_Matrix::M_Matrix(int num_states,int num_beads,double beta_num_beads)
    :num_states(num_states), num_beads(num_beads), beta_num_beads(beta_num_beads),

     M_alpha(num_states,num_states,0.0), V(num_states,num_beads,beta_num_beads),
     V_mat(num_beads,num_states,0.0),
     M_vec(num_beads,identity_matrix<std::complex<double> >(num_states)),
     exp_V_mat(num_beads,num_states,0.0),

     V_couple_vec(num_beads,identity_matrix<double>(num_states)),
     exp_row(num_states,0.0),V_couple_temp(num_states,num_states,0.0),
     M_inter(identity_matrix<double>(num_states))
{
    myExp.set_beta_num_beads(beta_num_beads);
}
void M_Matrix::update_M_vec(const vector<double> &Q){
        
    noalias(V_mat) = V.get_V_mat(Q);
    
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            exp_V_mat(bead,state) = exp(-beta_num_beads*V_mat(bead,state));
        }
    }
    
    update_V_couple_vec(Q);
    
    for (int bead=0; bead<num_beads; bead++) {
        update_M_alpha(Q[bead], bead);
        noalias(M_vec(bead)) = M_alpha;
    }
}
void M_Matrix::update_M_alpha(const double &Q, int bead){
    
    noalias(V_couple_temp) = V_couple_vec(bead);
    noalias(exp_row) = row(exp_V_mat, bead);
   
    for (int state1=0; state1<num_states-1; state1++) {
        for (int state2=1 + state1; state2<num_states; state2++) {
            M_inter(state1,state2) = -beta_num_beads*V_couple_temp(state1,state2);
            M_inter(state2,state1) = -beta_num_beads*V_couple_temp(state2,state1);
        }
    }

    for (int state=0; state<num_states; state++) {
        noalias(row(M_alpha,state)) = row(M_inter,state)*exp_row(state);
    }
}

void M_Matrix::exp_mat::set_beta_num_beads(double x){beta_num_beads = x;}

double M_Matrix::exp_mat::operator() (double x) const { return exp( - beta_num_beads * x); }

matrix<std::complex<double> >& M_Matrix::get_M_alpha(int alpha){
    return M_vec(alpha);
}
const matrix<double>& M_Matrix::get_V_mat() {return V_mat;}

void M_Matrix::update_V_couple_vec(const vector<double> &Q){
    
    for (int bead=0; bead<num_beads; bead++) {
        V_couple_vec(bead) = V.get_V_couple_mat(Q(bead));
    }
}
const matrix<double>& M_Matrix::get_V_couple_alpha(int alpha){
    return  V_couple_vec(alpha);
}

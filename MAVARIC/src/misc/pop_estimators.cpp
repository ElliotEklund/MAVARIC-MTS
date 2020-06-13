#include "pop_estimators.hpp"

pop_estimators::pop_estimators(int elec_beads,int num_states,double alpha)
    :elec_beads(elec_beads), num_states(num_states),
     x_sq(elec_beads,num_states,0),p_sq(elec_beads,num_states,0),
     xp_sq(elec_beads,num_states,0),
     x_sq_sum_bead(num_states,0),p_sq_sum_bead(num_states,0),
     x_sq_sum_state(elec_beads,0),p_sq_sum_state(elec_beads,0),
     one_elec_beads(elec_beads,1.0), one_num_states(num_states,1.0),
     one_mat(elec_beads,num_states,1.0),
     exp_v(elec_beads,0.0),
     wig_pops(num_states,0.0),
     gamma_nn(num_states,0.0),
     alpha(alpha)
{
    sc_coef = 1.0/(2.0*elec_beads);
    wigner_coef = pow(2.0,num_states+1)/elec_beads;
    x_alpha = alpha;
    p_alpha = 1.0/alpha;
}

vector<double> pop_estimators::boltz(const matrix<std::complex<double> > gamma){
    
    std::complex<double> tr(0.0);
    for (int state=0; state<num_states; state++) {
        tr += gamma(state,state);
        gamma_nn(state) = gamma(state,state).real();
    }
    
    return gamma_nn/tr.real();
}

vector<double> pop_estimators::sc(const matrix<double> &x, const matrix<double> &p){
    
    noalias(x_sq) = element_prod(x,x);
    noalias(p_sq) = element_prod(p,p);

    noalias(x_sq_sum_bead) = prod(one_elec_beads,x_sq);
    noalias(p_sq_sum_bead) = prod(one_elec_beads,p_sq);

    return sc_coef * (x_alpha*x_sq_sum_bead + p_alpha*p_sq_sum_bead -
                      elec_beads*one_num_states);
}

vector<double> pop_estimators::wigner(const matrix<double> &x, const matrix<double> &p){
    
    noalias(x_sq) = element_prod(x,x);
    noalias(p_sq) = element_prod(p,p);
    noalias(xp_sq) = x_sq + p_sq - 0.5*one_mat;

    noalias(x_sq_sum_state) = prod(x_sq,one_num_states);
    noalias(p_sq_sum_state) = prod(p_sq,one_num_states);
    
    for (int bead=0; bead<elec_beads; bead++) {
        exp_v(bead) = exp(-(x_sq_sum_state(bead) + p_sq_sum_state(bead)));
    }
    
    noalias(wig_pops) = prod(exp_v,xp_sq);
    
    return wigner_coef * wig_pops;
}

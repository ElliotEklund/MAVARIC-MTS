#include "dCdelec.hpp"

dCdelec::dCdelec(int elec_beads, int num_states, double alpha)
    :elec_beads(elec_beads), num_states(num_states),
     alpha(alpha),

     x_p_ip(elec_beads,num_states,0.0), x_m_ip(elec_beads,num_states,0.0),
     p_p_ix(elec_beads,num_states,0.0), p_m_ix(elec_beads,num_states,0.0),

     dCdx_row(num_states,0.0), dCdx_col(num_states,0.0),
     dCdp_row(num_states,0.0), dCdp_col(num_states,0.0),

     dCdx_mat(elec_beads,num_states,zero_matrix<std::complex<double> > (num_states,num_states)),
     dCdp_mat(elec_beads,num_states,zero_matrix<std::complex<double> > (num_states,num_states)),

     x_state(num_states,0.0), p_state(num_states,0.0),
     unit_complex(0.0,1.0)
{
    x_alpha = alpha;
    p_alpha = 1.0/alpha;
}
void dCdelec::update_dCdx_mat(const matrix<std::complex<double> > & x,
                              const matrix<std::complex<double> > & p){
    
//    noalias(x_p_ip) = x + unit_complex * p;
//    noalias(x_m_ip) = x - unit_complex * p;
    
    noalias(x_p_ip) = x_alpha*x + p_alpha * unit_complex * p;
    noalias(x_m_ip) = x_alpha*x - p_alpha * unit_complex * p;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            update_dCdx(bead, state, x,p, dCdx_mat(bead,state));
        }
    }    
}
void dCdelec::update_dCdp_mat(const matrix<std::complex<double> > & x,
                              const matrix<std::complex<double> > & p){
        
//    noalias(p_p_ix) = p + unit_complex * x;
//    noalias(p_m_ix) = p - unit_complex * x;
    
    noalias(p_p_ix) = p_alpha*p + x_alpha * unit_complex *x;
    noalias(p_m_ix) = p_alpha*p - x_alpha * unit_complex *x;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            update_dCdp(bead, state,x,p, dCdp_mat(bead,state));
        }
    }
}
void dCdelec::update_dCdx(int bead, int state, const matrix<std::complex<double> > & x,
                          const matrix<std::complex<double> > & p,
                          matrix<std::complex<double> > &dC){
    
    x_state = row(x,bead);
    p_state = row(p,bead);
    
    dCdx_row = x_alpha*x_state - unit_complex*p_state;
    dCdx_col = x_alpha*x_state + unit_complex*p_state;
    
    row(dC, state) = dCdx_row;
    column(dC, state) = dCdx_col;
    dC(state,state) = 2.0*x_alpha*x(bead,state);

//    dCdx_row = row(x_m_ip,bead);
//    dCdx_col = row(x_p_ip,bead);
//
//    row(dC, state) = dCdx_row;
//    column(dC, state) = dCdx_col;
//    dC(state,state) = 2.0*x;
}
void dCdelec::update_dCdp(int bead, int state, const matrix<std::complex<double> > & x,
                          const matrix<std::complex<double> > & p,
                          matrix<std::complex<double> > &dC){
    
    x_state = row(x,bead);
    p_state = row(p,bead);
    
    dCdp_row = p_alpha*p_state + unit_complex*x_state;
    dCdp_col = p_alpha*p_state - unit_complex*x_state;
    
    row(dC, state) = dCdp_row;
    column(dC, state) = dCdp_col;
    dC(state,state) = 2.0*p_alpha*p(bead,state);

//    dCdp_row = row(p_p_ix,bead);
//    dCdp_col = row(p_m_ix,bead);
//
//    row(dC, state) = dCdp_row;
//    column(dC, state) = dCdp_col;
//    dC(state,state) = 2.0*p;
}
const matrix<std::complex<double> > & dCdelec::get_dC_dx(int bead, int state){
    return dCdx_mat(bead,state);
}
const matrix<std::complex<double> > & dCdelec::get_dC_dp(int bead, int state){
    return dCdp_mat(bead,state);
}

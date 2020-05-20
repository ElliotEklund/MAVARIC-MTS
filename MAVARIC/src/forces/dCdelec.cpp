#include "dCdelec.hpp"

dCdelec::dCdelec(int elec_beads, int num_states)
    :elec_beads(elec_beads), num_states(num_states),

     x_p_ip(elec_beads,num_states,0.0), x_m_ip(elec_beads,num_states,0.0),
     p_p_ix(elec_beads,num_states,0.0), p_m_ix(elec_beads,num_states,0.0),

     dCdx_row(num_states,0.0), dCdx_col(num_states,0.0),
     dCdp_row(num_states,0.0), dCdp_col(num_states,0.0),

     dCdx_mat(elec_beads,num_states,zero_matrix<std::complex<double> > (num_states,num_states)),
     dCdp_mat(elec_beads,num_states,zero_matrix<std::complex<double> > (num_states,num_states)),

     unit_complex(0.0,1.0)

{}

void dCdelec::update_dCdx_mat(const matrix<std::complex<double> > & x,const matrix<std::complex<double> > & p){
    
    noalias(x_p_ip) = x + unit_complex * p;
    noalias(x_m_ip) = x - unit_complex * p;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            update_dCdx(bead, state, x(bead,state), dCdx_mat(bead,state));
        }
    }    
}

void dCdelec::update_dCdp_mat(const matrix<std::complex<double> > & x,const matrix<std::complex<double> > & p){
    
    noalias(p_p_ix) = p + unit_complex * x;
    noalias(p_m_ix) = p - unit_complex * x;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            update_dCdp(bead, state, p(bead,state), dCdp_mat(bead,state));
        }
    }
}

void dCdelec::update_dCdx(int bead, int state, std::complex<double> x, matrix<std::complex<double> > &dC){

    dCdx_row = row(x_m_ip,bead);
    dCdx_col = row(x_p_ip,bead);

    row(dC, state) = dCdx_row;
    column(dC, state) = dCdx_col;
    dC(state,state) = 2.0*x;
}

void dCdelec::update_dCdp(int bead, int state, std::complex<double> p, matrix<std::complex<double> > &dC){

    dCdp_row = row(p_p_ix,bead);
    dCdp_col = row(p_m_ix,bead);

    row(dC, state) = dCdp_row;
    column(dC, state) = dCdp_col;
    dC(state,state) = 2.0*p;
}

const matrix<std::complex<double> > & dCdelec::get_dC_dx(int bead, int state){
    return dCdx_mat(bead,state);
}

const matrix<std::complex<double> > & dCdelec::get_dC_dp(int bead, int state){
    return dCdp_mat(bead,state);
}

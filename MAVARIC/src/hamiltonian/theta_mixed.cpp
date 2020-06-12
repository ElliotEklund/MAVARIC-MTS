#include "theta_mixed.hpp"

theta_mixed::theta_mixed(int num_states,int nuc_beads, int elec_beads,
                         C_Matrix &C_In, M_Matrix &M_In)
    :num_states(num_states),
     nuc_beads(nuc_beads),
     elec_beads(elec_beads),
     gamma_mat(identity_matrix<std::complex<double> > (num_states)),
     Q_trans(elec_beads,0),
     W(elec_beads,nuc_beads,elec_beads)
{
    C = &C_In;
    M = &M_In;
    int ratio = nuc_beads/elec_beads;
    
    if (ratio == 0) {
        std::cout << "ERROR: nuc_beads is not divisible by elec_beads." << std::endl;
    }
    
    for (int i=0; i<elec_beads; i++) {
        W(i,i*ratio) = 1.0;
    }
}
void theta_mixed::update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                                   const matrix<double> &p){
    
    Q_trans = prod(W,Q);
    C->update_C_vec(x, p);
    M->update_M_vec(Q_trans);
        
    /* Chain multiplicaton */
    gamma_mat = identity_matrix<std::complex<double> > (num_states);

    for(int bead=0; bead<elec_beads; bead++){
        gamma_mat = prod(gamma_mat,C->get_C_alpha(bead));
        gamma_mat = prod(gamma_mat,M->get_M_alpha(bead));
    }
}
void theta_mixed::update_theta(const vector<double> &Q,const matrix<double> &x,
                               const matrix<double> &p){
    
    update_gamma_mat(Q, x, p);
    std::complex<double> tr (0.0,0.0);
    tr = trace<std::complex<double> >(gamma_mat,num_states);
    theta = tr.real();
}
double theta_mixed::get_theta(const vector<double> &Q,const matrix<double> &x,
                              const matrix<double> &p){
    update_theta(Q, x, p);
    return theta;
}
double theta_mixed::get_theta(){return theta;}

double theta_mixed::get_signTheta(){return sign(theta);}

double theta_mixed::get_signTheta(const vector<double> &Q,const matrix<double> &x,
                                  const matrix<double> &p){
    
    update_theta(Q,x,p);
    return sign(theta);
}
matrix<std::complex<double> > theta_mixed::get_gamm(){return gamma_mat;}

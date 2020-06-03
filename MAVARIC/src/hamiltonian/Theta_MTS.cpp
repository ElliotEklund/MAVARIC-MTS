#include "Theta_MTS.hpp"

Theta_MTS::Theta_MTS(int num_states, int elec_beads,
                     C_Matrix &C_In, M_Matrix_MTS &M_MTS_In)
    :num_states(num_states),
     elec_beads(elec_beads),
     gamma_mat(identity_matrix<std::complex<double> > (num_states))

{
    C = &C_In;
    M_MTS = &M_MTS_In;
}
void Theta_MTS::update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                             const matrix<double> &p){
    
    C->update_C_vec(x, p);
    M_MTS->update_M_MTS_vec(Q);
        
    /* Chain multiplicaton */
    gamma_mat = identity_matrix<std::complex<double> > (num_states);

    for(int bead=0; bead<elec_beads; bead++){
        gamma_mat = prod(gamma_mat,C->get_C_alpha(bead));
        gamma_mat = prod(gamma_mat,M_MTS->get_M_MTS_alpha(bead));
    }
}
void Theta_MTS::update_theta(const vector<double> &Q,const matrix<double> &x,
                         const matrix<double> &p){
    
    update_gamma_mat(Q, x, p);
    std::complex<double> tr (0.0,0.0);
    tr = trace<std::complex<double> >(gamma_mat,num_states);
    theta = tr.real();
}
double Theta_MTS::get_theta(const vector<double> &Q,const matrix<double> &x,
                        const matrix<double> &p){
    
    update_theta(Q, x, p);
    return theta;
}
double Theta_MTS::get_theta(){return theta;}

double Theta_MTS::get_signTheta(){return sign(theta);}

double Theta_MTS::get_signTheta(const vector<double> &Q,const matrix<double> &x,
                                const matrix<double> &p){
    
    update_theta(Q,x,p);
    return sign(theta);
}
matrix<std::complex<double> > Theta_MTS::get_gamm(){return gamma_mat;}

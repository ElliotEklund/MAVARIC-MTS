#include "mvrpmd_mixed_forces.hpp"

mvrpmd_mixed_forces::mvrpmd_mixed_forces(int nuc_beads,int elec_beads, int num_states,
                                         double mass, double beta_nuc_beads, double alpha,
                                         theta_mixed &thetaIN, theta_mixed_dQ &theta_dQIN,
                                         theta_mixed_dElec &theta_dElec_IN)

    :ONE_mass(1.0/mass),beta_nuc_beads(beta_nuc_beads),
     nuc_beads(nuc_beads), ONE_beta_nuc_beads(1.0/beta_nuc_beads),
     TWO_beta_nuc_beads(2.0/beta_nuc_beads),elec_beads(elec_beads),
     num_states(num_states),

     dVspring_dQ(nuc_beads, mass, beta_nuc_beads),
     dV0_dQ(nuc_beads,mass),

     dHdQ(nuc_beads,0.0), dHdP(nuc_beads,0.0),
     dHdx(elec_beads,num_states,0.0),dHdp(elec_beads,num_states,0.0),

     dVspring_dQ_vec(nuc_beads,0.0), dV0_dQ_vec(nuc_beads,0.0),
     dThetaMTS_dQ_vec(nuc_beads,0.0),

     dThetaMTS_dx_vec(elec_beads,num_states,0.0),
     dThetaMTS_dp_vec(elec_beads,num_states,0.0),

     alpha(alpha)
{
    theta = &thetaIN;
    theta_dQ = &theta_dQIN;
    theta_dElec = &theta_dElec_IN;
    coeff_ONE_theta = 0;
    
    x_alpha = alpha;
    p_alpha = 1.0/alpha;
}
void mvrpmd_mixed_forces::update_Forces(const vector<double> &Q,const vector<double> &P,
                                        const matrix<double> &x,const matrix<double> &p){

    /* Calling update_theta updates C and M_MTS. All further function calls that
     used C and M_MTS do not need to update these objects. */
    theta->update_theta(Q, x, p);
    theta_dQ->update_theta_dQ(Q);
    theta_dElec->update_theta_dElec(x, p);

    coeff_ONE_theta = ONE_beta_nuc_beads/theta->get_theta();

    update_dHdP(P);
    update_dHdQ(Q, x, p);
    update_dHdx(Q, x, p);
    update_dHdp(Q, x, p);
}
void mvrpmd_mixed_forces::update_dHdP(const vector<double> &P){
    dHdP =  ONE_mass * P;
}
void mvrpmd_mixed_forces::update_dHdQ(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){
    
    dVspring_dQ_vec = dVspring_dQ.get_dSpring_dQ(Q);
    dV0_dQ_vec = dV0_dQ.get_dStateIndep_dQ(Q);
    dThetaMTS_dQ_vec = theta_dQ->get_theta_dQ_vec();
    noalias(dHdQ) = dVspring_dQ_vec + dV0_dQ_vec - coeff_ONE_theta*dThetaMTS_dQ_vec;
}
void mvrpmd_mixed_forces::update_dHdx(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){
    dThetaMTS_dx_vec = theta_dElec->get_theta_dx_vec();
    //noalias(dHdx) = x_alpha*TWO_beta_nuc_beads*x - coeff_ONE_theta * dThetaMTS_dx_vec;

    noalias(dHdx) = /*TWO_beta_nuc_beads*x*/ - coeff_ONE_theta * dThetaMTS_dx_vec;
}
void mvrpmd_mixed_forces::update_dHdp(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){

    dThetaMTS_dp_vec = theta_dElec->get_theta_dp_vec();
    //noalias(dHdp) = p_alpha*TWO_beta_nuc_beads*p - coeff_ONE_theta * dThetaMTS_dp_vec;

    noalias(dHdp) = /*TWO_beta_nuc_beads*p*/ - coeff_ONE_theta * dThetaMTS_dp_vec;
}
const vector<double> & mvrpmd_mixed_forces::get_dHdQ(){return dHdQ;}

const vector<double> & mvrpmd_mixed_forces::get_dHdP(){return dHdP;}

const matrix<double> & mvrpmd_mixed_forces::get_dHdx(){return dHdx;}

const matrix<double> & mvrpmd_mixed_forces::get_dHdp(){return dHdp;}

double mvrpmd_mixed_forces::get_sign(const vector<double> &Q, const matrix<double> &x,
                                const matrix<double> &p){
    return theta->get_signTheta(Q,x,p);
}
void mvrpmd_mixed_forces::print_dHdQ(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdQ" << std::endl;
    std::cout << dHdQ << std::endl;
    std::cout << std::endl << std::endl;
}
void mvrpmd_mixed_forces::print_dHdP(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdP" << std::endl;
    std::cout << dHdP << std::endl;
    std::cout << std::endl << std::endl;
}
void mvrpmd_mixed_forces::print_dHdx(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdx" << std::endl;
    std::cout << dHdx << std::endl;
    std::cout << std::endl << std::endl;
}
void mvrpmd_mixed_forces::print_dHdp(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdp" << std::endl;
    std::cout << dHdp << std::endl;
    std::cout << std::endl << std::endl;
}

#include "Forces_MTS.hpp"

Forces_MTS::Forces_MTS(int nuc_beads,int elec_beads, int num_states,
                       double mass, double beta_nuc_beads,
                       Theta_MTS &ThetaMTS_In, dTheta_MTS_dQ &dThetaMTS_dQ_In,
                       dTheta_MTS_dElec &dThetaMTS_dElec_In)

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
     dThetaMTS_dp_vec(elec_beads,num_states,0.0)

{
    ThetaMTS = &ThetaMTS_In;
    dThetaMTS_dQ = &dThetaMTS_dQ_In;
    dThetaMTS_dElec = &dThetaMTS_dElec_In;

    coeff_ONE_theta = 0;
}
void Forces_MTS::update_Forces(const vector<double> &Q,const vector<double> &P,
                   const matrix<double> &x,const matrix<double> &p){

    /* Calling update_theta updates C and M_MTS. All further function calls that
     used C and M_MTS do not need to update these objects. */
    ThetaMTS->update_theta(Q, x, p);
    dThetaMTS_dQ->update_dTheta_MTS_dQ_vec(Q);
    dThetaMTS_dElec->update_dTheta_MTS_dElec(x, p);

    coeff_ONE_theta = ONE_beta_nuc_beads/ThetaMTS->get_theta();

    update_dHdP(P);
    update_dHdQ(Q, x, p);
    update_dHdx(Q, x, p);
    update_dHdp(Q, x, p);
}

void Forces_MTS::update_dHdP(const vector<double> &P){
    dHdP =  ONE_mass * P;
    //dHdP = nuc_beads* ONE_mass * P;
}

void Forces_MTS::update_dHdQ(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){
    
    dVspring_dQ_vec = dVspring_dQ.get_dSpring_dQ(Q);
    dV0_dQ_vec = dV0_dQ.get_dStateIndep_dQ(Q);
    dThetaMTS_dQ_vec = dThetaMTS_dQ->get_dThetaMTS_dQ_vec();
    noalias(dHdQ) = dVspring_dQ_vec + dV0_dQ_vec - coeff_ONE_theta*dThetaMTS_dQ_vec;
    
    /* N pulled out verision */
    //noalias(dHdQ) = (dVspring_dQ_vec + dV0_dQ_vec - coeff_ONE_theta*dThetaMTS_dQ_vec)/nuc_beads;
}

void Forces_MTS::update_dHdx(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){
    
    dThetaMTS_dx_vec = dThetaMTS_dElec->get_dThetaMTS_dx_vec();
    noalias(dHdx) = /*(TWO_beta_nuc_beads*x*/ - coeff_ONE_theta * dThetaMTS_dx_vec;
    //noalias(dHdx) = - coeff_ONE_theta * dThetaMTS_dx_vec/nuc_beads;
}

void Forces_MTS::update_dHdp(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){

    dThetaMTS_dp_vec = dThetaMTS_dElec->get_dThetaMTS_dp_vec();
    noalias(dHdp) = /*(TWO_beta_nuc_beads*p*/ - coeff_ONE_theta * dThetaMTS_dp_vec;
    //noalias(dHdp) = - coeff_ONE_theta * dThetaMTS_dp_vec/nuc_beads;
}

const vector<double> & Forces_MTS::get_dHdQ(){return dHdQ;}

const vector<double> & Forces_MTS::get_dHdP(){return dHdP;}

const matrix<double> & Forces_MTS::get_dHdx(){return dHdx;}

const matrix<double> & Forces_MTS::get_dHdp(){return dHdp;}

double Forces_MTS::get_sign(const vector<double> &Q, const matrix<double> &x,
                                const matrix<double> &p){
    return ThetaMTS->get_signTheta(Q,x,p);
}
void Forces_MTS::print_dHdQ(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdQ" << std::endl;
    std::cout << dHdQ << std::endl;
    std::cout << std::endl << std::endl;
}

void Forces_MTS::print_dHdP(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdP" << std::endl;
    std::cout << dHdP << std::endl;
    std::cout << std::endl << std::endl;
}

void Forces_MTS::print_dHdx(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdx" << std::endl;
    std::cout << dHdx << std::endl;
    std::cout << std::endl << std::endl;
}

void Forces_MTS::print_dHdp(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdp" << std::endl;
    std::cout << dHdp << std::endl;
    std::cout << std::endl << std::endl;
}

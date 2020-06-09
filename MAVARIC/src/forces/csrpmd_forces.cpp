#include "csrpmd_forces.hpp"

csrpmd_forces::csrpmd_forces(int nuc_beads,int elec_beads, int num_states,
                                         double mass, double beta_nuc_beads)
    :ONE_mass(1.0/mass),
     beta_nuc_beads(beta_nuc_beads),
     nuc_beads(nuc_beads),
     elec_beads(elec_beads),
     num_states(num_states),

     dVspring_dQ(nuc_beads, mass, beta_nuc_beads),
     dV0_dQ(nuc_beads,mass),
     Vsc(nuc_beads,elec_beads,num_states),

     dVspring_dQ_vec(nuc_beads,0.0),
     dV0_dQ_vec(nuc_beads,0.0),
     Vsc_dQ_vec(nuc_beads,0.0),

     dHdQ(nuc_beads,0.0), dHdP(nuc_beads,0.0),
     dHdx(elec_beads,num_states,0.0),dHdp(elec_beads,num_states,0.0)
{
    
}
void csrpmd_forces::update_Forces(const vector<double> &Q,const vector<double> &P,
                   const matrix<double> &x,const matrix<double> &p){
    update_dHdP(P);
    update_dHdQ(Q, x, p);
    update_dHdx(Q,x);
    update_dHdp(Q,p);
}
void csrpmd_forces::update_dHdP(const vector<double> &P){
    noalias(dHdP) =  ONE_mass * P;
}
void csrpmd_forces::update_dHdQ(const vector<double> &Q, const matrix<double> &x,
                                       const matrix<double> &p){

    noalias(dVspring_dQ_vec) = dVspring_dQ.get_dSpring_dQ(Q);
    noalias(dV0_dQ_vec) = dV0_dQ.get_dStateIndep_dQ(Q);
    noalias(Vsc_dQ_vec) = Vsc.get_Vsc_dQ(Q,x,p);
    noalias(dHdQ) = dVspring_dQ_vec + dV0_dQ_vec + Vsc_dQ_vec;
}
void csrpmd_forces::update_dHdx(const vector<double> &Q, const matrix<double> &x){
    noalias(dHdx) = Vsc.get_Vsc_dx(Q,x);
}
void csrpmd_forces::update_dHdp(const vector<double> &Q, const matrix<double> &p){
    noalias(dHdp) = Vsc.get_Vsc_dp(Q,p);
}
double csrpmd_forces::get_sign(const vector<double> &Q, const matrix<double> &x,
                               const matrix<double> &p){
    return 1.0;
}
const vector<double> & csrpmd_forces::get_dHdQ(){return dHdQ;}

const vector<double> & csrpmd_forces::get_dHdP(){return dHdP;}

const matrix<double> & csrpmd_forces::get_dHdx(){return dHdx;}

const matrix<double> & csrpmd_forces::get_dHdp(){return dHdp;}

void csrpmd_forces::print_dHdQ(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdQ" << std::endl;
    std::cout << dHdQ << std::endl;
    std::cout << std::endl << std::endl;
}
void csrpmd_forces::print_dHdP(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdP" << std::endl;
    std::cout << dHdP << std::endl;
    std::cout << std::endl << std::endl;
}
void csrpmd_forces::print_dHdx(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdx" << std::endl;
    std::cout << dHdx << std::endl;
    std::cout << std::endl << std::endl;
}
void csrpmd_forces::print_dHdp(){
    std::cout << std::endl << std::endl;
    std::cout << "dHdp" << std::endl;
    std::cout << dHdp << std::endl;
    std::cout << std::endl << std::endl;
}

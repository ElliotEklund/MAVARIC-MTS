#include "Forces_two_particles.hpp"

Forces_two_particles::Forces_two_particles(int nuc_beads1, int nuc_beads2, double mass, double beta)

    :nuc_beads1(nuc_beads1),nuc_beads2(nuc_beads2),
     W(nuc_beads1,nuc_beads2),

     dVspring_dQ1(nuc_beads1,mass,beta/nuc_beads1),
     dVspring_dQ2(nuc_beads2,mass,beta/nuc_beads2),

     dHdQ1(nuc_beads1,0.0),dHdQ2(nuc_beads2,0.0),
     dHdP1(nuc_beads1,0.0),dHdP2(nuc_beads2,0.0),
 
     dVspring_dQ1_vec(nuc_beads1,0.0),
     dVspring_dQ2_vec(nuc_beads2,0.0),

     dV0_dQ1_vec(nuc_beads1,0.0),
     dV0_dQ2_vec(nuc_beads2,0.0),
     dV12_dQ1_vec(nuc_beads1,0.0),
     dV12_dQ2_vec(nuc_beads2,0.0)

{
    c1 = nuc_beads1*1.0/(mass);
    c2 = nuc_beads2*1.0/(mass);
    
    d1 = 0.5*(5.0/nuc_beads1);
    d2 = 0.5*(1.0/nuc_beads2);
    d12 = 0*(2.0/nuc_beads1);

    one_nuc_beads1 = 1.0/nuc_beads1;
    one_nuc_beads2 = 1.0/nuc_beads2;
    
    int r = nuc_beads1/nuc_beads2;
    
    if (nuc_beads1==1 && nuc_beads2==1) {
        W(0,0) = 0.0;
    }
    else{
        for (int i=0; i<nuc_beads2; i++) {
            for (int j=0; j<r; j++) {
                W(i*r+j,i) = 1.0;
            }
        }
    }
}

void Forces_two_particles::update_Forces(const vector<double> &Q1,const vector<double> &Q2,
                                         const vector<double> &P1,const vector<double> &P2){
    update_dHdP1(P1);
    update_dHdP2(P2);
    update_dHdQ1(Q1,Q2);
    update_dHdQ2(Q1,Q2);
}
void Forces_two_particles::update_dHdP1(const vector<double> &P1){
    noalias(dHdP1) = c1 * P1;
}
void Forces_two_particles::update_dHdP2(const vector<double> &P2){
    noalias(dHdP2) = c2 * P2;
}
void Forces_two_particles::update_dHdQ1(const vector<double> &Q1, const vector<double> &Q2){
    noalias(dVspring_dQ1_vec) = one_nuc_beads1*dVspring_dQ1.get_dSpring_dQ2(Q1);
    noalias(dV0_dQ1_vec) = one_nuc_beads1*dV0_dQ1(Q1);
    noalias(dV12_dQ1_vec) = dV12_dQ1(Q1,Q2);
    noalias(dHdQ1) = dVspring_dQ1_vec + dV0_dQ1_vec + dV12_dQ1_vec;
}
void Forces_two_particles::update_dHdQ2(const vector<double> &Q1, const vector<double> &Q2){
    noalias(dVspring_dQ2_vec) = one_nuc_beads2*dVspring_dQ2.get_dSpring_dQ2(Q2);
    noalias(dV0_dQ2_vec) = one_nuc_beads2*dV0_dQ2(Q2);
    noalias(dV12_dQ2_vec) = dV12_dQ2(Q1,Q2);
    noalias(dHdQ2) = dVspring_dQ2_vec + dV0_dQ2_vec + dV12_dQ2_vec;
}
vector<double> Forces_two_particles::dV0_dQ1(const vector<double> &Q){

    vector<double> QQ = element_prod(Q,Q);
    vector<double> QQQ = element_prod(QQ,Q);
    
    return Q + 0.3*QQ + 0.04*QQQ;
}
vector<double> Forces_two_particles::dV0_dQ2(const vector<double> &Q){
    vector<double> QQ = element_prod(Q,Q);
    vector<double> QQQ = element_prod(QQ,Q);
   
    return  QQQ;
}
vector<double> Forces_two_particles::dV12_dQ1(const vector<double> &Q1,const vector<double> &Q2){
    return 2.0*d12*(Q1-prod(W,Q2));
}
vector<double> Forces_two_particles::dV12_dQ2(const vector<double> &Q1,const vector<double> &Q2){
    return 2.0*d12* prod((prod(W,Q2) - Q1),W);
}
const vector<double> & Forces_two_particles::get_dHdQ1(){return dHdQ1;}

const vector<double> & Forces_two_particles::get_dHdQ2(){return dHdQ2;}

const vector<double> & Forces_two_particles::get_dHdP1(){return dHdP1;}

const vector<double> & Forces_two_particles::get_dHdP2(){return dHdP2;}

void Forces_two_particles::print_dHdQ(int sub_system){
    std::cout << std::endl << std::endl;
    std::cout << "dHdQ" << sub_system << std::endl;
    if (sub_system==1) {
        std::cout << dHdQ1 << std::endl;
    }
    else {
        std::cout << dHdQ2 << std::endl;
    }
    std::cout << std::endl << std::endl;
}
void Forces_two_particles::print_dHdP(int sub_system){
    std::cout << std::endl << std::endl;
    std::cout << "dHdP" << sub_system << std::endl;
    if (sub_system==1) {
        std::cout << dHdP1 << std::endl;
    }
    else {
        std::cout << dHdP2 << std::endl;
    }
    std::cout << std::endl << std::endl;
}

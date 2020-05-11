#include "two_particle_Hamiltonian.hpp"

two_particle_Hamiltonian::two_particle_Hamiltonian(int num_beads1, int num_beads2,
                                                   SpringEnergy &V_springIn1,SpringEnergy &V_springIn2)
    :num_beads1(num_beads1),num_beads2(num_beads2),
     one_num_beads1(1.0/num_beads1),one_num_beads2(1.0/num_beads2),
     W(num_beads1,num_beads2),
     Q2_mapped(num_beads1,0.0)

{
    V_spring1 = &V_springIn1;
    V_spring2 = &V_springIn2;
    
    c1 = 5.0/num_beads1;
    c2 = 1.0/num_beads2;
    c12 = 0.0/num_beads1;
    
    int r = num_beads1/num_beads2;
    
    for (int i=0; i<num_beads2; i++) {
        for (int j=0; j<r; j++) {
            W(i*r+j,i) = 1.0;
        }
    }
}
double inline two_particle_Hamiltonian::V01(const vector<double> &Q1){
    vector<double> QQ = element_prod(Q1,Q1);
    vector<double> QQQ = element_prod(QQ,Q1);
    vector<double> QQQQ = element_prod(QQQ,Q1);

    double t_2 = sum(QQ);
    double t_3 = sum(QQQ);
    double t_4 = sum(QQQQ);

    return (0.5*t_2 + 0.1*t_3 + 0.01*t_4)/num_beads1;
}
double inline two_particle_Hamiltonian::V02(const vector<double> &Q2){
    
    vector<double> QQ = element_prod(Q2,Q2);
    vector<double> QQQQ = element_prod(QQ,QQ);
    double t_4 = sum(QQQQ);

    return 0.25*t_4/num_beads2;
}
double inline two_particle_Hamiltonian::Vcouple(const vector<double> &Q1,const vector<double> &Q2){
    noalias(Q2_mapped) = prod(W,Q2);
    return c12* inner_prod(Q1-Q2_mapped,Q1-Q2_mapped);
}
double two_particle_Hamiltonian::get_energy(const vector<double> &Q1, const vector<double> &Q2){
        
    double E_spring1 = one_num_beads1*V_spring1->get_springEnergy(Q1);
    double E_spring2 = one_num_beads2*V_spring2->get_springEnergy(Q2);
    double E_V01 = V01(Q1);
    double E_V02 = V02(Q2);
    double E_coup = Vcouple(Q1,Q2);
    return E_spring1 + E_spring2 + E_V01 + E_V02 + E_coup;
}
//
//double two_particle_Hamiltonian::get_energy(){
//
//    double V_spring_temp = V_spring->get_springEnergy();
//    double V0_temp = V0->get_V0();
//
//    return V_spring_temp + V0_temp;
//}

//double two_particle_Hamiltonian::get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P){
//
//    double kin_energy = 0.5 * inner_prod(P,P)/mass;
//    double V_spring_temp = V_spring->get_springEnergy(Q);
//    double V0_temp = V0->get_V0(Q);
//
//    return kin_energy + V_spring_temp + V0_temp;
//}

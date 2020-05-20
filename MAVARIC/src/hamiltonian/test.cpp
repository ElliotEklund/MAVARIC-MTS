#include <iostream>
#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "rpmd_vv.hpp"
#include "RK4_MVRPMD.hpp"
#include "Forces_MTS.hpp"
#include "mvrpmd_special.hpp"

using namespace boost::numeric::ublas;

void fill_P(vector<double> &P){
    P(0) = -1.18386;
    P(1) =  0.985709;
    P(2) = -3.01999;
    P(3) = -0.015654;
    P(4) = -2.76225;
    P(5) = -1.46936;
}
void fill_Q(vector<double> &Q){
    Q(0) = -0.373708;
    Q(1) =-0.31678;
    Q(2) =-0.356903;
    Q(3) =-0.521235;
    Q(4) =-0.686323;
    Q(5) =-0.717537;
}
void fill_xElec(matrix<double> &x){
    x(0,0) = -1.37131;
    x(0,1) =0.0461138;
    x(1,0) =1.05452;
    x(1,1) =0.688434;
    x(2,0) =0.77922;
    x(2,1) =-1.19172;
    x(3,0) =-1.99172;
    x(3,1) =0.930253;
    x(4,0) =-0.543918;
    x(4,1) =0.478252;
    x(5,0) =-0.259731;
    x(5,1) =-2.27023;
}
void fill_pElec(matrix<double> &p){
    p(0,0) = 1.45801;
    p(0,1) =-0.61309;
    p(1,0) =-0.475432;
    p(1,1) =-0.35186;
    p(2,0) =-1.10597;
    p(2,1) =-0.534775;
    p(3,0) =-1.7657;
    p(3,1) =-1.29253;
    p(4,0) =-0.70437;
    p(4,1) =-2.04044;
    p(5,0) =1.17777;
    p(5,1) =0.263807;
}
void write(const vector<double> &v, std::ofstream & myStream){

    int s = v.size();
    
    for (int i=0; i<s; i++) {
        myStream << v[i] << " ";
    }
    myStream << std::endl;
}

void write(const matrix<double> &m, std::ofstream & myStream){
    
    int num_row = m.size1();
    int num_col = m.size2();
    
    for (int row=0; row<num_row; row++) {
        for (int col=0; col<num_col; col++) {
            myStream << m(row,col) << " ";
        }
    }

    myStream << std::endl;
}

int main(int argc, char ** argv){

    int nuc_beads, elec_beads, num_states;
    double beta, mass;
    double dt, time;
    
    nuc_beads = 6;
    elec_beads = 6;
    num_states = 2;
    beta = 1.0;
    mass = 1.0;
    
    dt = 0.001;
    time = 10;
    
    vector<double> Q(nuc_beads,0), P(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0), p(elec_beads,num_states,0);
    
    fill_Q(Q);
    fill_P(P);
    fill_xElec(x);
    fill_pElec(p);
    
    P = P/nuc_beads;
    x = x/nuc_beads;
    p = p/nuc_beads;
    
    C_Matrix C(elec_beads, num_states);
    M_Matrix M(num_states, nuc_beads, beta/nuc_beads);
    M_Matrix_MTS M_MTS(nuc_beads, elec_beads, num_states, M);
    dM_Matrix_dQ dMdQ(nuc_beads, num_states, beta/elec_beads, M);
    dM_Matrix_MTS_dQ dM_MTS_dQ(nuc_beads, elec_beads, num_states, dMdQ);
    Theta_MTS Theta(num_states, elec_beads, C, M_MTS);
    dTheta_MTS_dQ dThetadQ(num_states, nuc_beads, elec_beads, C, M_MTS, dM_MTS_dQ);
    dTheta_MTS_dElec dThetadElec(num_states, elec_beads, C, M_MTS);
    
    Forces_MTS F(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads,
                 Theta,dThetadQ,dThetadElec);
    
    Forces_MTS * Fp = &F;
    
    RK4_MVRPMD rk4(Fp,nuc_beads,elec_beads,num_states,dt);
    
    int num_steps = time/dt;
    
    std::ofstream Qstream, Pstream, xStream, pStream;
    std::string root = "Data/";
    std::string tag = "H7";
    
    std::string Q_file, P_file, xElec_file, pElec_file;
    Q_file = root + "Q_" + tag;
    P_file = root + "P_" + tag;
    xElec_file = root + "xElec_" + tag;
    pElec_file = root + "pElec_" + tag;
    
    Qstream.open(Q_file);
    Pstream.open(P_file);
    xStream.open(xElec_file);
    pStream.open(pElec_file);

    for (int step=0; step<num_steps; step++) {
        rk4.take_step(Q,P,x,p);
        write(Q,Qstream);
        write(P,Pstream);
        write(x,xStream);
        write(p,pStream);
    }
    
    Qstream.close();
    Pstream.close();
    xStream.close();
    pStream.close();

    return 0;
}

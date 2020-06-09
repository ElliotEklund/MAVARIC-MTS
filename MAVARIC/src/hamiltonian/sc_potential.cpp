#include "sc_potential.hpp"

sc_potential::sc_potential(int nuc_beads, int elec_beads, int num_states)
    :nuc_beads(nuc_beads),
     elec_beads(elec_beads),
     num_states(num_states),
     W(elec_beads,nuc_beads,elec_beads),
     V(nuc_beads,elec_beads,elec_beads),

     Q_trans(elec_beads,0),
     x_temp(num_states,0),p_temp(num_states,0),
     Vsc_dQ_temp(elec_beads,0),
     Vsc_dQ(nuc_beads,0),
     Vsc_dx(elec_beads,num_states,0),Vsc_dp(elec_beads,num_states,0)
{
    pot = 0.0;
    double delta = 1.0; //constant coupling parameter
    Vmat.resize(num_states,num_states);
    Vmat_dQ.resize(num_states,num_states);

    for (int i=0; i<num_states; i++) {
        for (int j=0; j<num_states; j++) {
            Vmat(i,j) = delta;
            Vmat_dQ(i,j) = 0.0;
        }
    }

    int ratio = nuc_beads/elec_beads;
    
    if (ratio == 0) {
        std::cout << "ERROR: nuc_beads is not divisible by elec_beads." << std::endl;
    }
    
    for (int i=0; i<elec_beads; i++) {
        W(i,i*ratio) = 1.0;
        V(i*ratio,i) = 1.0;
    }
}
void sc_potential::update_Vsc(const vector<double> &Q,const matrix<double> &x,
                              const matrix<double> &p){

    Q_trans = prod(W,Q);
    double sum = 0;
    
    for (int bead=0; bead<elec_beads; bead++) {
        matrix_row <const matrix<double> > x_row (x, bead);
        matrix_row <const matrix<double> > p_row (p, bead);

        update_Vmat(Q_trans(bead),Vmat);
                
        noalias(x_temp) = prod(Vmat,x_row);
        noalias(p_temp) = prod(Vmat,p_row);

        sum += 0.5*(inner_prod(x_row,x_temp) + inner_prod(p_row,p_temp));
        sum -= trace(Vmat,num_states);
    }
    pot = sum;
}
void sc_potential::update_Vsc_dQ(const vector<double> &Q,const matrix<double> &x,
                                 const matrix<double> &p){

    Q_trans = prod(W,Q);
    
    for (int bead=0; bead<elec_beads; bead++) {
        matrix_row <const matrix<double> > x_row (x, bead);
        matrix_row <const matrix<double> > p_row (p, bead);

        update_Vmat_dQ(Q_trans(bead),Vmat_dQ);
                
        noalias(x_temp) = prod(Vmat_dQ,x_row);
        noalias(p_temp) = prod(Vmat_dQ,p_row);

        Vsc_dQ_temp(bead) = 0.5*(inner_prod(x_row,x_temp) + inner_prod(p_row,p_temp));
        Vsc_dQ_temp(bead) -= trace(Vmat_dQ,num_states);
    }
    Vsc_dQ = prod(V,Vsc_dQ_temp);
}
void sc_potential::update_Vsc_dx(const vector<double> &Q,const matrix<double> &x){
    
    Q_trans = prod(W,Q);
    for (int bead=0; bead<elec_beads; bead++) {
        
        matrix_row <matrix<double> > Vsc_dx_row (Vsc_dx, bead);
        matrix_row <const matrix<double> > x_row (x, bead);

        update_Vmat(Q_trans(bead),Vmat);        
        Vsc_dx_row = prod(Vmat,x_row);
    }
}
void sc_potential::update_Vsc_dp(const vector<double> &Q,const matrix<double> &p){
    
    Q_trans = prod(W,Q);
    for (int bead=0; bead<elec_beads; bead++) {
        matrix_row <matrix<double> > Vsc_dp_row (Vsc_dp, bead);
        matrix_row <const matrix<double> > p_row (p, bead);

        update_Vmat(Q_trans(bead),Vmat);
        Vsc_dp_row = prod(Vmat,p_row);
    }
}
inline double sc_potential::V_d(const double &Q, int n){
    if (n==0) {return Q;}
    else{return -Q;}
}
inline double sc_potential::V_d_dQ(const double &Q, int n){
    if (n==0) {return 1.0;}
    else{return -1.0;}
}
inline double sc_potential::V_od(const double &Q, int n, int m){
    return 1;
}
inline double sc_potential::V_od_dQ(const double &Q, int n, int m){
    return 0.0;
}
inline void sc_potential::update_Vmat(const double &Q, matrix<double> &VmatIN){
    for (int i=0; i<num_states; i++){
        VmatIN(i,i) = V_d(Q,i);
    }
}
inline void sc_potential::update_Vmat_dQ(const double &Q, matrix<double> &Vmat_dQIN){
    for (int i=0; i<num_states; i++){
        Vmat_dQIN(i,i) = V_d_dQ(Q,i);
    }
}
double sc_potential::get_Vsc(){return pot;}

double sc_potential::get_Vsc(const vector<double> &Q,const matrix<double> &x,
                             const matrix<double> &p){
    update_Vsc(Q,x,p);
    return pot;
}
const matrix<double> & sc_potential::get_Vsc_dx(const vector<double> &Q,
                                                const matrix<double> &p){
    
    update_Vsc_dx(Q,p);
    return Vsc_dx;
}
const matrix<double> & sc_potential::get_Vsc_dp(const vector<double> &Q,
                                                const matrix<double> &x){
    
    update_Vsc_dp(Q,x);
    return Vsc_dp;
}

const vector<double> & sc_potential::get_Vsc_dQ(const vector<double> &Q,
                                                const matrix<double> &x,
                                                const matrix<double> &p){
    
    update_Vsc_dQ(Q,x,p);
    return Vsc_dQ;
}

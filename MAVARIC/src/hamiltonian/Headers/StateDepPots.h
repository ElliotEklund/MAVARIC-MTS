#ifndef STATEDEPPOTS_H
#define STATEDEPPOTS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <algorithm>

using namespace boost::numeric::ublas;

class StateDepPots{
  
public:
    StateDepPots(int num_states, int num_beads, double beta_num_beads);
    
    /* Return a vector where each element is (V11(Q1),...V11(QN))*/
    vector<double> get_V11_vec(const vector<double> &Q);
   
    /* Return a vector where each element is (V22(Q1),...V22(QN))*/
    vector<double> get_V22_vec(const vector<double> &Q);
    
    vector<double> get_V33_vec(const vector<double> &Q);

    /* Return a matrix where the diagonal is 1; all other elements
     are Vi,j(Q)*/
    matrix<double>& get_V_couple_mat(const double &Q);

    /* Return a matrix where Mij = Vjj(Qi)*/
    matrix<double>& get_V_mat(const vector<double> &Q);
    
private:
    
    /* Private Structs. */
    
    /* Apply V11 to all elements */
    struct V11 {
        double operator() (double x) const;
    };
    
    /* Apply V22 to all elements */
    struct V22 {
        double operator() (double x) const;
    };

    /* Apply V33 to all elements */
    struct V33 {
        double operator() (double x) const;
    };
    
    double V12(const double &Q);

    double V13(const double &Q);
  
    double V23(const double &Q);

    /* Private Data*/
    int num_states; //number of electronic states
    int num_beads; //number of ring polymer beads
    double beta_num_beads; //beta/num_beads
    
    matrix<double> V_mat; //Vii matrix evaluated over all Q
    matrix<double> V_couple_mat; //matrix of V12 evaluated over all Q

    vector<double> V11_vec;
    vector<double> V22_vec;
    vector<double> V33_vec;

};

#endif

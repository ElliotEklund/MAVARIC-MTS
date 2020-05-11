#ifndef VV_TWO_PARTICLE_HPP
#define VV_TWO_PARTICLE_HPP

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "Forces_two_particles.hpp"
#include "dSpring_dQ.hpp"

using namespace boost::numeric::ublas;


class VV_two_particle{
    
public:
    
    VV_two_particle(int nuc_beads1, int nuc_beads2, double mass, double beta,
                    double dt);
    
    void step(vector<double> &Q1,vector<double> &Q2,
              vector<double> &P1,vector<double> &P2);

private:
    
    double dt;
    double HALF_dt;
    vector<double> P1_half, P2_half;
    
    Forces_two_particles F;
    
};

#endif

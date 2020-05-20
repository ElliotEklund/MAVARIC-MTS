#ifndef DECORRELATION_H
#define DECORRELATION_H

#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

using namespace boost::numeric::ublas;


class decorrelation{
    
public:
    decorrelation(int nuc_beads, int elec_beads, int num_samples);
    
    /* Return nuclear centroid given Q*/
    double centroid(const vector<double> &Q);
    
    /* Compute the centroid of Q and save it in all_centroids[i] */
    void stash_samples(int i, const vector<double> &Q);

    /* Write all_centroids to file*/
    void write_all_centroids(std::string root);
    
private:
    int nuc_beads; //number of beads,
    double one_nuc_beads; //1.0/nuc_beads
    int elec_beads; //number of beads,
    int num_samples; //number of samples used
    vector<double> one; //vector of ones
    vector<double> all_centroids; //collection of centroids

    
    
};


#endif

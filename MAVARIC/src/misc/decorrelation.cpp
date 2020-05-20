#include "decorrelation.h"

decorrelation::decorrelation(int nuc_beads, int elec_beads, int num_samples)
    :nuc_beads(nuc_beads), one_nuc_beads(1.0/nuc_beads),
     elec_beads(elec_beads),
     num_samples(num_samples),
     one(nuc_beads,1.0),
     all_centroids(num_samples,0.0)
{}

double decorrelation::centroid(const vector<double> &Q){
    return one_nuc_beads * inner_prod(one,Q);
}

void decorrelation::stash_samples(int i, const vector<double> &Q){
    all_centroids(i) = centroid(Q);
}

void decorrelation::write_all_centroids(std::string root){
    
    std::string fileName = root + "Output/decorrelation";
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if (!myFile.is_open()){
        std::cout << "ERROR: Could not open file " << fileName << std::endl;
    }
    
    for (int i=0; i<num_samples; i++) {
        myFile << all_centroids[i] << std::endl;
    }
    
    myFile.close();
}

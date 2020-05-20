#ifndef SAMPLING_H
#define SAMPLING_H

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>

#include "errors.h"

/*
 Return 1 if no errors were found in Sampling inputfile.
 Else, return 0.
 */
inline int check_sampling_in(std::string fileName, std::vector<double> &params){
    
    int result = 0;
    std::ifstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    std::string myLine;
    
    while(getline(myFile,myLine)){
        int foundAt = myLine.find(":"); //location of colon
        std::string goodStuff = myLine.substr(foundAt + 1, myLine.size());
        params.push_back(atof(goodStuff.c_str()));
    }
    
    myFile.close();
    
    /* Check that Run Sampling is binary. */
    if(!is_binary(params[0],fileName,"Run Sampling")){
        result = -1;
    }
    
    /* Check that Number of Trajectories is a positive integer. */
    if(!is_positive_int(params[1],fileName,"Number of Trajectories")){
        result = -1;
    }
    
    /* Check that Decorrelation Length is a positive integer. */
    if(!is_positive_int(params[2],fileName,"Decorrelation Length")){
        result = -1;
    }
    
    /* Check that Save Sampled Trajectories is binary. */
    if(!is_binary(params[3],fileName,"Save Sampled Trajectories")){
        result = -1;
    }
    /* Check that Histogram Positions is binary. */
    if(!is_binary(params[4],fileName,"Histogram Positions")){
        result = -1;
    }
    
    /* Check that Number of Bins is a positive integer. */
    if(!is_positive_int(params[5],fileName,"Number of Bins")){
        result = -1;
    }
    
    /* Check that Read PSV is binary. */
    if(!is_binary(params[6],fileName,"Read PSV")){
        result = -1;
    }
    
    return result;
}

#endif

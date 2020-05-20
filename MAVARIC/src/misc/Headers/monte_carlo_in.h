#ifndef MONTE_CARLO_IN_H
#define MONTE_CARLO_IN_H

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>

#include "errors.h"

/*
 Return 1 if no errors were found in MonteCarlo inputfile.
 Else, return 0.
 */
inline int check_monte_carlo_in(std::string fileName, std::vector<double> &params){
    
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
    
    /* Check that Run MC is binary. */
    if(!is_binary(params[0],fileName,"Run MC")){
        result = -1;
    }
    
    /* Check that MC Moves ia positive integer. */
    if(!is_positive_int(params[1],fileName,"MC Moves")){
        result = -1;
    }
    
    /* Check that Estimator Rate ia positive integer. */
    if(!is_positive_int(params[2],fileName,"Estimator Rate")){
        result = -1;
    }
    
    /* Check that Estimator Rate is smaller than MC Moves */
    if(params[2] >= params[1]){
        std::cout << "ERROR: In file MonteCarlo, Estimator Rate must be less than MC Moves."
        " Aborting calculation." << std::endl;
        return -1;
    }
    
    /* Check that MC Moves is divisible by Estimator Rate. */
    if(int(params[1]) % int(params[2]) != 0){
        std::cout << "ERROR: In file MonteCarlo, MC Moves must be divisible by Estimator Rate."
        " Aborting calculation." << std::endl;
        return -1;
    }
    
    /*Check that Save PSV is binary */
    if(!is_binary(params[3],fileName,"Save PSV")){
        result = -1;
    }
    
    /*Check that Save MC Data is binary */
    if(!is_binary(params[4],fileName,"Save MC Data")){
        result = -1;
    }
    
    /*Check that Read PSV is binary */
    if(!is_binary(params[5],fileName,"Read PSV")){
        result = -1;
    }
    
    /*Check that Read MC Data is binary */
    if(!is_binary(params[6],fileName,"Read PSV")){
        result = -1;
    }
    
    return result;
}

#endif

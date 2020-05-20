#ifndef DYNAMICS_IN_H
#define DYNAMICS_IN_H

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>

#include "errors.h"

/*
 Return 1 if no errors were found in Dynamics inputfile.
 Else, return 0.
 */
inline int check_dynamics_in(std::string fileName, std::vector<double> &params){
    
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
    
    /* Check that Run Dynamics is binary. */
    if(!is_binary(params[0],fileName,"Run Dynamics")){
        result = -1;
    }
    
    /* Check that Run Time is a positive number. */
    if(!is_positive(params[1],fileName,"Run Time")){
        result = -1;
    }
    
    /* Check that Time Step is a positive number. */
    if(!is_positive(params[2],fileName,"Time Step")){
        result = -1;
    }
    
    double intpart; //used to check for integer
    
    /* Check that Run Time is divisible by Time Step */
    if(!(modf(params[1]/params[2],&intpart) == 0.0)){
        std::cout << "ERROR: In file Dynamics, Run Time must be divisible by Time Step."
        " Aborting calculation." << std::endl;
        return -1;
    }
    
    /* Check that Read Trajectories is binary. */
    if(!is_binary(params[3],fileName,"Read Trajectories")){
        result = -1;
    }
    
    /* Check that Compute CQQ is binary. */
    if(!is_binary(params[4],fileName,"Compute CQQ")){
        result = -1;
    }
    
    /* Check that Check Energy Conservation is binary. */
    if(!is_binary(params[5],fileName,"Check Energy Conservation")){
        result = -1;
    }
    
    /* Check that Conservation Tolerance is a positive number. */
    if(!is_positive(params[6],fileName,"Conservation Tolerance")){
        result = -1;
    }
    
    /* Check that Compute Wigner Populations is binary. */
    if(!is_binary(params[7],fileName,"Compute Wigner Populations")){
        result = -1;
    }

    /* Check that Compute Initial PAC is binary. */
    if(!is_binary(params[8],fileName,"Compute Initial PAC")){
        result = -1;
    }

    /* Check that Interval is a positive number. */
    if(!is_positive(params[9],fileName,"Interval")){
        result = -1;
    }
    
    return result;
}

#endif

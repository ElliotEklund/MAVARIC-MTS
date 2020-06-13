#ifndef ELECTRONIC_IN_H
#define ELECTRONIC_IN_H

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>

#include "errors.h"

/*
 Return 1 if no errors were found in ElectronicParameters inputfile.
 Else, return 0.
 */
inline int check_elec_in(std::string root, std::vector<double> &params){
    
    int result = 0;
    
    std::ifstream myFile;
    std::string fileName = root;
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
    
    /* Check Number of States */
    if(!is_positive_int(params[0],fileName,"States")){
        result = -1;
    }
    
    if(!is_positive_int(params[1],fileName,"Beads")){
        result = -1;
    }
    
    /* Check x step Size */
    if(!is_positive(params[2],fileName,"x step size")){
        result = -1;
    };
    
    /* Check p step Size */
    if(!is_positive(params[3],fileName,"p step size")){
        result = -1;
    };
    
    /* Check alpha */
    if(!is_positive(params[4],fileName,"alpha")){
        result = -1;
    };
    
    /* If no errors have been caught, return 0. */
    return result;
}
#endif

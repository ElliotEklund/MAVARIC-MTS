#ifndef SYSTEM_IN
#define SYSTEM_IN

#include "errors.h"

#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>

/*
 Return 1 if no errors were found in SystemParameters inputfile.
 Else, return 0.
 */
inline int check_sys_in(std::string fileName, std::vector<double> & params){
    
    int result = 0;
    std::ifstream myFile;
    myFile.open(fileName.c_str());

    if(!myFile.is_open()) {
        std::cout << "ERROR: Could not open file " << fileName << std::endl;
    }

    std::string myLine;
    
    while(getline(myFile,myLine)){
        int foundAt = myLine.find(":"); //location of colon
        std::string goodStuff = myLine.substr(foundAt + 1, myLine.size());
        params.push_back(atof(goodStuff.c_str()));
    }

    myFile.close();

    /* Check Mass */
    if(!is_positive(params[0],fileName,"Mass")){result = -1;}

    /* Check Beads */
    if(!is_positive_int(params[1],fileName,"Beads")){result = -1;}

    /* Check Temperature */
    if(!is_positive(params[2],fileName,"Temperature")){result = -1;}

    /* Check MC Step Size */
    if(!is_positive(params[3],fileName,"MC Step Size")){result = -1;}

     //If no errors have been caught, return 0.
    return result;
}

#endif

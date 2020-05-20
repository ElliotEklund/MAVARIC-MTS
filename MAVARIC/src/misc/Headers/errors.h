#ifndef ERRORS_H
#define ERRORS_H

#include <iostream>
#include <string>
#include <math.h>
#include <fstream>

/*
 Return true if x is a positive number, else return false.
 If false, an error message will be printed.
*/
inline bool is_positive(double x, std::string fileName, std::string varName){
    
    if (x < 0){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive number."
        " Aborting calculation." << std::endl;
        return false;
    }
    else{
        return true;
    }
}

/*
 Return true if x is a positive int, else return false.
 If false, an error message will be printed.
 */
inline bool is_positive_int(double x, std::string fileName, std::string varName){

    if ((x < 0) || (x != floor(x))){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive integer."
        " Aborting calculation." << std::endl;
        return false;
    }
    else{
        return true;
    }
}

/*
Return true if x is binary (0,1), else return false.
If false, an error message will be printed.
*/
inline bool is_binary(double x, std::string fileName, std::string varName){

    if((x!=0) && (x!=1)){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be either 1 or 0."
        " Aborting calculation." << std::endl;
        return false;
    }
    else{
        return true;
    }
}

/*
 Return true if the file is not corrupt, else return false.
 A file is corrupt if it is missing a colon, does not contain a
 phrase before the colon, or is missing a input after the colon.
 */
inline bool is_corrupt(std::string fileName, int correct_num_lines ){

    std::ifstream myFile;
    myFile.open(fileName.c_str());

    if(!myFile.is_open()) {
        std::cout << "ERROR: Could not open file " << fileName << std::endl;
    }

    std::string myLine;
    int num_lines = 0;
    bool file_corrupt = false;

    while(getline(myFile, myLine)){

        num_lines += 1;

        if(!myLine.empty()){

            int foundAt = myLine.find(":"); //location of colon
            std::string goodStuff = myLine.substr(foundAt + 1, myLine.size());

            if((foundAt <= 0) || (goodStuff.length()==0)){
                file_corrupt = true;
            }

            if(goodStuff.find(' ') != -1){
                file_corrupt = true;
            }

            char* p;
            double converted = strtod(goodStuff.c_str(), &p);

            if (*p) {
                file_corrupt = true;
            }
        }
    }

    if(num_lines != correct_num_lines){
        file_corrupt = true;
    }

    myFile.close();

    if(file_corrupt){
        std::cout << "ERROR: File " << fileName << " is corrupt. Aborting calculation." << std::endl;
    }

    return file_corrupt;
}

#endif

#ifndef MainHlpr_hpp
#define MainHlpr_hpp

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

class MainHlpr{
    
    
public:
    MainHlpr();
    
    /* Read in the parameters from the InputFile "fileName". Return
     * these parameters in a vector.*/
    std::vector<double> get_parameters(std::string fileName);

    /* corruption_check returns true if fileName is corrupt.
     * A file is corrupt if it meets any of the following criteria.
     * 1. The number of lines in fileName does not match correct_num_lines.
     * 2. The format DESCRIPTION:VALUE is altered.
     * 3. Any white space or non-numerical characters are found after the colon.*/
    bool corruption_check(std::string fileName, int correct_num_lines );

    /* Return true if x is a positive number. Otherwise, return false. */
    bool check_positive(double x, std::string fileName, std::string varName);

    /* Return true if x is a positive integer. Otherwise, return false. */
    bool check_positive_int(double x, std::string fileName, std::string varName);

    /* Return true if x i either 0 or 1. Otherwise, return false. */
    bool check_binary(double x, std::string fileName, std::string varName);

    /* Calling this function assumes that SysParameters is not corrupt and passed
     * the call to corruption_check. */
    int SysParams_checker(std::string root, std::vector<double> & params);
    
    /* Calling this function assumes that ElecParameters is not corrupt and passed
     * the call to corruption_check. */
    int ElecParams_checker(std::string root, std::vector<double> &params);

    /* Calling this function assumes that MonteCarlo is not corrupt and passed
     * the call to corruption check. */
    int MCParams_checker(std::string root, std::vector<double> &params);

    /* Calling this function assumes that Sampling is not corrupt and passed
     * the call to corruption check. */
    int SampParams_checker(std::string root, std::vector<double> &params);

    /* Calling this function assumes that Dynamics is not corrupt and passed
     * the call to corruption check. */
    int DynParams_checker(std::string root, std::vector<double> &params);

    /* Each of the five input vectors well be filled with their corresponding values after
     * input_file_handler is called. These vectors will only be returned if the function does
     * not find any corrupt files or errors.*/
    int input_file_handler(std::string root, std::vector<double> &sysParams,std::vector<double> &elecParams,
                                  std::vector<double> &MCParams,std::vector<double> &sampParams,
                                  std::vector<double> &dynParams);
    
};

#endif /* MainHlpr_hpp */

#include "input_rpmd_bath.hpp"

input_rpmd_bath::input_rpmd_bath()
{}

std::vector<double> input_rpmd_bath::get_parameters(std::string fileName){

    std::vector<double> params;

    std::ifstream myFile;
    myFile.open(fileName.c_str());

    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }

    std::string myLine;

    while(getline(myFile, myLine)){

        if(!myLine.empty()){
            size_t foundAt = myLine.find(":"); //location of colon
            std::string goodStuff = myLine.substr(foundAt + 1, myLine.size());
            params.push_back(atof(goodStuff.c_str()));
        }
        else { /* no statement */ }
    }

    myFile.close();
    return params;
}


/* Each of the five input vectors well be filled with their corresponding values after
 * input_file_handler is called. These vectors will only be returned if the function does
 * not find any corrupt files or errors.*/
int input_rpmd_bath::input_file_handler(std::string root, std::vector<double> &sysParams,
                                        std::vector<double> &elecParams,std::vector<double> &MCParams,
                                        std::vector<double> &sampParams,std::vector<double> &dynParams){
    
    
    std::string ElecFile = root + "InputFiles/ElecParameters.txt";
    std::string SysFile = root + "InputFiles/SystemParameters.txt";
    std::string MCFile = root + "InputFiles/MonteCarlo.txt";
    std::string SampFile = root + "InputFiles/Sampling.txt";
    std::string DynFile = root + "InputFiles/Dynamics.txt";
    
    int result = 0;

    /* First, check if any of the input files are corrupt */
    if(is_corrupt(ElecFile,3)){
        result = -1;
    }

    if(is_corrupt(SysFile,4)){
        result = -1;
    }

    if(is_corrupt(MCFile,7)){
        result = -1;
    }

    if(is_corrupt(SampFile,7)){
        result = -1;
    }

    if(is_corrupt(DynFile,10)){
        result = -1;
    }

    /* If an input file is corrupt, do not proceed exit
     * input_file_handler and return -1 */
    if(result == -1){
        return -1;
    }

    /* Check that the parameters in each file are acceptable. */
    if(check_sys_in(root, sysParams) == -1){
        result = -1;
    };

    if(check_elec_in(root, elecParams) == -1){
        result = -1;
    };

    if(check_monte_carlo_in(root, MCParams) == -1){
        result = -1;
    }

    if(check_sampling_in(root, sampParams) == -1){
        result = -1;
    }

    if(check_dynamics_in(root, dynParams) == -1){
        result = -1;
    }
    
    return result;
}

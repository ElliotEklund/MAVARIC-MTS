#include "MainHlpr.hpp"

MainHlpr::MainHlpr()
{}

std::vector<double> MainHlpr::get_parameters(std::string fileName){

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

bool MainHlpr::corruption_check(std::string fileName, int correct_num_lines ){

    std::ifstream myFile;
    myFile.open(fileName.c_str());

    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
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

bool MainHlpr::check_positive(double x, std::string fileName, std::string varName){

    if (x < 0){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive number."
        " Aborting calculation." << std::endl;
        return false;
    }
    else{
        return true;
    }
}

bool MainHlpr::check_positive_int(double x, std::string fileName, std::string varName){

    if ((x < 0) || (x != floor(x))){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive integer."
        " Aborting calculation." << std::endl;
        return false;
    }
    else{
        return true;
    }
}

bool MainHlpr::check_binary(double x, std::string fileName, std::string varName){

    if((x!=0) && (x!=1)){
        std::cout << "ERROR: In file " + fileName + ", " + varName + " must be either 1 or 0."
        " Aborting calculation." << std::endl;
        return false;
    }

    else{
        return true;
    }
}

int MainHlpr::SysParams_checker(std::string root, std::vector<double> & params){
    
    int result = 0;
    std::ifstream myFile;
    std::string fileName = root + "InputFiles/SystemParameters.txt";
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
    
    /* Check Mass */
    if(!check_positive(params[0],fileName,"Mass")){
        result = -1;
    }
    
    /* Check Beads */
    if(!check_positive_int(params[1],fileName,"Beads")){
        result = -1;
    }
    
    /* Check Temperature */
    if(!check_positive(params[2],fileName,"Temperature")){
        result = -1;
    }
    
    /* Check MC Step Size */
    if(!check_positive(params[3],fileName,"MC Step Size")){
        result = -1;
    }
    
    /* If no errors have been caught, return 0. */
    return result;
}

/* Calling this function assumes that ElecParameters is not corrupt and passed
 * the call to corruption_check. */
int MainHlpr::ElecParams_checker(std::string root, std::vector<double> &params){
    
    int result = 0;
    
    std::ifstream myFile;
    std::string fileName = root + "InputFiles/ElecParameters.txt";
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
    if(!check_positive_int(params[0],fileName,"States")){
        result = -1;
    }
    
    if(!check_positive_int(params[1],fileName,"Beads")){
        result = -1;
    }
    
    /* Check MC Step Size */
    if(!check_positive(params[2],fileName,"MC Step Size")){
        result = -1;
    };
    
    /* If no errors have been caught, return 0. */
    return result;
}

/* Calling this function assumes that MonteCarlo is not corrupt and passed
 * the call to corruption check. */
int MainHlpr::MCParams_checker(std::string root, std::vector<double> &params){
    
    int result = 0;
    std::ifstream myFile;
    std::string fileName = root + "InputFiles/MonteCarlo.txt";
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
    if(!check_binary(params[0],fileName,"Run MC")){
        result = -1;
    }
    
    /* Check that MC Moves ia positive integer. */
    if(!check_positive_int(params[1],fileName,"MC Moves")){
        result = -1;
    }
    
    /* Check that Estimator Rate ia positive integer. */
    if(!check_positive_int(params[2],fileName,"Estimator Rate")){
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
    if(!check_binary(params[3],fileName,"Save PSV")){
        result = -1;
    }
    
    /*Check that Save MC Data is binary */
    if(!check_binary(params[4],fileName,"Save MC Data")){
        result = -1;
    }
    
    /*Check that Read PSV is binary */
    if(!check_binary(params[5],fileName,"Read PSV")){
        result = -1;
    }
    
    /*Check that Read MC Data is binary */
    if(!check_binary(params[6],fileName,"Read PSV")){
        result = -1;
    }
    
    return result;
}

/* Calling this function assumes that Sampling is not corrupt and passed
 * the call to corruption check. */
int MainHlpr::SampParams_checker(std::string root, std::vector<double> &params){
    
    int result = 0;
    std::ifstream myFile;
    std::string fileName = root + "InputFiles/Sampling.txt";
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
    if(!check_binary(params[0],fileName,"Run Sampling")){
        result = -1;
    }
    
    /* Check that Number of Trajectories is a positive integer. */
    if(!check_positive_int(params[1],fileName,"Number of Trajectories")){
        result = -1;
    }
    
    /* Check that Decorrelation Length is a positive integer. */
    if(!check_positive_int(params[2],fileName,"Decorrelation Length")){
        result = -1;
    }
    
    /* Check that Save Sampled Trajectories is binary. */
    if(!check_binary(params[3],fileName,"Save Sampled Trajectories")){
        result = -1;
    }
    /* Check that Histogram Positions is binary. */
    if(!check_binary(params[4],fileName,"Histogram Positions")){
        result = -1;
    }
    
    /* Check that Number of Bins is a positive integer. */
    if(!check_positive_int(params[5],fileName,"Number of Bins")){
        result = -1;
    }
    
    /* Check that Read PSV is binary. */
    if(!check_binary(params[6],fileName,"Read PSV")){
        result = -1;
    }
    
    return result;
}

/* Calling this function assumes that Dynamics is not corrupt and passed
 * the call to corruption check. */
int MainHlpr::DynParams_checker(std::string root, std::vector<double> &params){
    
    int result = 0;
    std::ifstream myFile;
    std::string fileName = root + "InputFiles/Dynamics.txt";
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
    if(!check_binary(params[0],fileName,"Run Dynamics")){
        result = -1;
    }
    
    /* Check that Run Time is a positive number. */
    if(!check_positive(params[1],fileName,"Run Time")){
        result = -1;
    }
    
    /* Check that Time Step is a positive number. */
    if(!check_positive(params[2],fileName,"Time Step")){
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
    if(!check_binary(params[3],fileName,"Read Trajectories")){
        result = -1;
    }
    
    /* Check that Compute CQQ is binary. */
    if(!check_binary(params[4],fileName,"Compute CQQ")){
        result = -1;
    }
    
    /* Check that Check Energy Conservation is binary. */
    if(!check_binary(params[5],fileName,"Check Energy Conservation")){
        result = -1;
    }
    
    /* Check that Conservation Tolerance is a positive number. */
    if(!check_positive(params[6],fileName,"Conservation Tolerance")){
        result = -1;
    }
    
    /* Check that Compute Wigner Populations is binary. */
    if(!check_binary(params[7],fileName,"Compute Wigner Populations")){
        result = -1;
    }

    /* Check that Compute Initial PAC is binary. */
    if(!check_binary(params[8],fileName,"Compute Initial PAC")){
        result = -1;
    }

    /* Check that Interval is a positive number. */
    if(!check_positive(params[9],fileName,"Interval")){
        result = -1;
    }
    
    return result;
}

/* Each of the five input vectors well be filled with their corresponding values after
 * input_file_handler is called. These vectors will only be returned if the function does
 * not find any corrupt files or errors.*/
int MainHlpr::input_file_handler(std::string root, std::vector<double> &sysParams,std::vector<double> &elecParams,
                              std::vector<double> &MCParams,std::vector<double> &sampParams,
                              std::vector<double> &dynParams){
    
    
    std::string ElecFile = root + "InputFiles/ElecParameters.txt";
    std::string SysFile = root + "InputFiles/SystemParameters.txt";
    std::string MCFile = root + "InputFiles/MonteCarlo.txt";
    std::string SampFile = root + "InputFiles/Sampling.txt";
    std::string DynFile = root + "InputFiles/Dynamics.txt";
    
    int result = 0;
    
    /* First, check if any of the input files are corrupt */
    if(corruption_check(ElecFile,3)){
        result = -1;
    }
    
    if(corruption_check(SysFile,4)){
        result = -1;
    }
    
    if(corruption_check(MCFile,7)){
        result = -1;
    }
    
    if(corruption_check(SampFile,7)){
        result = -1;
    }
    
    if(corruption_check(DynFile,10)){
        result = -1;
    }
    
    /* If an input file is corrupt, do not proceed exit
     * input_file_handler and return -1 */
    if(result == -1){
        return -1;
    }
    
    
    /* Check that the parameters in each file are acceptable. */
    if(SysParams_checker(root, sysParams) == -1){
        result = -1;
    };
    
    if(ElecParams_checker(root, elecParams) == -1){
        result = -1;
    };
    
    if(MCParams_checker(root, MCParams) == -1){
        result = -1;
    }
    
    if(SampParams_checker(root, sampParams) == -1){
        result = -1;
    }
    
    if(DynParams_checker(root, dynParams) == -1){
        result = -1;
    }
    
    return result;
}

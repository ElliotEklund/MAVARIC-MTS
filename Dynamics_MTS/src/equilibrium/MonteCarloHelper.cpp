#include "MonteCarloHelper.h"

MonteCarloHelper::MonteCarloHelper(std::string root)
    :root(root)
{}

void MonteCarloHelper::print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt){
    std::cout << "\t System Acceptance Ratio: " << 100*(double)sys_steps_accpt/sys_steps << std::endl;
}

void MonteCarloHelper::print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt){
    std::cout << "\t Electronic Acceptance Ratio: " << 100*(double)elec_steps_accpt/elec_steps << std::endl;
}

void MonteCarloHelper::print_avg_energy(double estimator_total, double sgn_total){
    std::cout << "\t Average Energy: " << estimator_total/sgn_total << std::endl;
}

void MonteCarloHelper::write_estimator(vector <double> estimator,int interval){
    
    int size = estimator.size();
    
    std::string fileName = root + "/Results/estimator.txt";
    
    std::ofstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    for(int i=0; i<size; i++){
        myFile << i*interval << " " <<  estimator[i] << std::endl;
    }
    
    myFile.close();
    
    std::cout << "\t Successfully wrote energy_estimator file to Results." << std::endl;
}

void MonteCarloHelper::write_PSV(int nuc_beads,int elec_beads, int num_states,
                      vector<double> Q, matrix<double> x,matrix <double> p){
    
    std::string fileName = root + "/Results/PSV.txt";
    
    std::ofstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    for(int bead=0; bead<nuc_beads; bead++){
        myFile << Q(bead) << std::endl;
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile << x(bead,state) << std::endl;
        }
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile << p(bead,state) << std::endl;
        }
    }
    
    std::cout << "\t Successfully saved PSV to Results." << std::endl;
    myFile.close();
}

void MonteCarloHelper::read_PSV(int nuc_beads,int elec_beads, int num_states, vector<double> &Q, matrix<double> &x,
                     matrix<double> &p){
    
    std::string fileName = root + "/Results/PSV.txt";
    
    std::ifstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    for(int bead=0; bead<nuc_beads; bead++){
        myFile >> Q(bead);
    }
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile >> x(bead,state);
        }
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile >> p(bead,state);
        }
    }
    
    myFile.close();
    
    std::cout << "\t Successfully read PSV file from Results." << std::endl;
}

void MonteCarloHelper::write_MC_data(double sgn_total,double estimator_total){
    
    std::string fileName = root + "/Results/mcData.txt";
    
    std::ofstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile << sgn_total << std::endl;
    myFile << estimator_total << std::endl;
    
    myFile.close();
    
    std::cout << "\t Successfully saved MC data to Results." << std::endl;
    
}

void MonteCarloHelper::read_MC_data(double &sgn_totalGlobal, double &estimator_total){
    
    std::string fileName = root + "/Results/mcData.txt";
    
    std::ifstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile >> sgn_totalGlobal;
    myFile >> estimator_total;
    
    myFile.close();
    
    std::cout << "\t Successfully read MC datat from Results." << std::endl;
}

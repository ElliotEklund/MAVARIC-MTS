#include "MonteCarloHelper.h"

MonteCarloHelper::MonteCarloHelper(std::string root, int my_id, int num_procs, int root_proc)
    :root(root), my_id(my_id),
     num_procs(num_procs), root_proc(root_proc)
{
    sys_ratio = 0;
    x_ratio = 0;
    p_ratio = 0;
    average_energy = 0;
    energy_stdev = 0;
}
void MonteCarloHelper::set_sys_ratio(unsigned long long sys_steps,
                                     unsigned long long sys_steps_accpt){
    
    unsigned long long sys_steps_global = 0;
    unsigned long long sys_steps_accpt_global = 0;
    
    MPI_Reduce(&sys_steps,&sys_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
               root_proc, MPI_COMM_WORLD);
    
    MPI_Reduce(&sys_steps_accpt,&sys_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    sys_ratio =  100*(double)sys_steps_accpt_global/sys_steps_global;
}
void MonteCarloHelper::set_x_ratio(unsigned long long x_steps,
                                   unsigned long long x_steps_accpt){
    
    unsigned long long x_steps_global = 0;
    unsigned long long x_steps_accpt_global = 0;
    
    MPI_Reduce(&x_steps,&x_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
               root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&x_steps_accpt,&x_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    x_ratio = 100*(double)x_steps_accpt_global/x_steps_global;
}
void MonteCarloHelper::set_p_ratio(unsigned long long p_steps,
                                   unsigned long long p_steps_accpt){
    
    unsigned long long p_steps_global = 0;
    unsigned long long p_steps_accpt_global = 0;
    
    MPI_Reduce(&p_steps,&p_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
               root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&p_steps_accpt,&p_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    p_ratio = 100*(double)p_steps_accpt_global/p_steps_global;
}
void MonteCarloHelper::print_sys_accpt(){
    
    if (my_id == root_proc) {
        std::cout << "\t System Acceptance Ratio: " << sys_ratio << std::endl;
    }
}
void MonteCarloHelper::print_x_accpt(){
    
    if (my_id == root_proc) {
        std::cout << "\t x Acceptance Ratio: " << x_ratio<< std::endl;
    }
}
void MonteCarloHelper::print_p_accpt(){
    
    if (my_id == root_proc) {
        std::cout << "\t p Acceptance Ratio: " << p_ratio<< std::endl;
    }
}
void MonteCarloHelper::set_average_energy(double estimator_total, double sgn_total){
    
    double estimator_total_global = 0;
    double sgn_total_global = 0;
    double standev = 0;
    double standev_global = 0;
    
    MPI_Reduce(&estimator_total,&estimator_total_global,1,MPI_DOUBLE,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&sgn_total,&sgn_total_global,1,MPI_DOUBLE,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        estimator_total_global = estimator_total_global/sgn_total_global;
    }
    
    MPI_Bcast(&estimator_total_global,1,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    standev = pow(estimator_total_global - (estimator_total/sgn_total),2);
    MPI_Reduce(&standev,&standev_global,1,MPI_DOUBLE,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    average_energy = estimator_total_global;
    
    if (num_procs > 1) {
        energy_stdev = sqrt(standev_global/(num_procs-1));}
    else{
        energy_stdev = 0;}
}
void MonteCarloHelper::print_avg_energy(){

    if (my_id == root_proc) {
        if (num_procs > 1) {
            std::cout << "\t Average Energy: " << average_energy << std::endl;
            std::cout << "\t Standard Deviation: " << energy_stdev << std::endl;
            std::cout << "\t Using " << num_procs << " processors." << std::endl;
            std::cout << "\t STDEV/Estimator:" << 100*energy_stdev/average_energy
            << "%" << std::endl;
        }
        else{
            std::cout << "\t Average Energy: " << average_energy << std::endl;
            std::cout << "\t Standard Deviation: NA " << std::endl;
        }
    }
}
void MonteCarloHelper::write_estimator(vector <double> estimator,int interval){
    
    int size = estimator.size();
    
    std::string fileName = root + "Output/estimator.txt";
    
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    for(int i=0; i<size; i++){
        myFile << i*interval << " " <<  estimator[i] << std::endl;
    }
    
    myFile.close();
    
    if (my_id==root_proc) {
        std::cout << "\t Successfully wrote energy_estimator file to Results." << std::endl;
    }
}
void MonteCarloHelper::write_PSV(int nuc_beads,int elec_beads, int num_states,
                      vector<double> Q, matrix<double> x,matrix <double> p){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "Output/PSV" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
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
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully saved PSV to Results." << std::endl;
    }
    
    myFile.close();
}
void MonteCarloHelper::read_PSV(int nuc_beads,int elec_beads, int num_states, vector<double> &Q, matrix<double> &x,
                     matrix<double> &p){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "Output/PSV" + quickConvert.str();
    
    std::ifstream myFile;
    myFile.open(fileName.c_str());
    
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
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully read PSV file from Results." << std::endl;
    }
}
void MonteCarloHelper::write_MC_data(double sgn_total,double estimator_total){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "Output/mcData" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile << sgn_total << std::endl;
    myFile << estimator_total << std::endl;
    
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully saved MC data to Results." << std::endl;
    }
    
}
void MonteCarloHelper::read_MC_data(double &sgn_totalGlobal, double &estimator_total){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "Output/mcData" + quickConvert.str();
    
    std::ifstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile >> sgn_totalGlobal;
    myFile >> estimator_total;
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully read MC datat from Results." << std::endl;
    }
}
void MonteCarloHelper::final_report(int nuc_beads, int elec_beads, int num_states,
                                    double beta, unsigned long long num_steps,
                                    double nuc_ss, double x_ss, double p_ss){
    
    if (my_id == root_proc) {
        std::string file_name = root + "Output/equil_report";

        std::ofstream myStream;
        myStream.open(file_name.c_str());
        
        if(!myStream.is_open()) {
            std::cout << "ERROR: Could not open file" << file_name << std::endl;
        }
    
        myStream << "nuc_beads:" << nuc_beads << std::endl;
        myStream << "elec_beads:" << elec_beads << std::endl;
        myStream << "num_states:" << num_states << std::endl;
        myStream << "beta:" << beta << std::endl;
        myStream << "Monte Carlo Steps:" << num_steps << std::endl;
        myStream << "Energy Estimator:" << average_energy << std::endl;
        
        if (num_procs > 1) {
            myStream << "Energy Estimator STDEV:" << energy_stdev << std::endl;
            myStream << "STDEV/Estimator:" << 100* energy_stdev/average_energy <<
            "%" << std::endl;
        }
        else{
            myStream << "Energy Estimator STDEV: NA" << std::endl;
        }

        myStream << "Nuclear Step Size:" << nuc_ss << " Acceptance Ratio: " <<
        sys_ratio << std::endl;
        
        myStream << "x Step Size:" << x_ss << " Acceptance Ratio: " <<
        x_ratio << std::endl;
        
        myStream << "p Step Size:" << p_ss << " Acceptance Ratio: " <<
        p_ratio << std::endl;
        
        myStream.close();
    }
    
    print_sys_accpt();
    print_x_accpt();
    print_p_accpt();
    print_avg_energy();
}

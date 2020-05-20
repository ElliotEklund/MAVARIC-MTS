#include "MonteCarloHelper.h"

MonteCarloHelper::MonteCarloHelper(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs), root_proc(root_proc)
{}

void MonteCarloHelper::set_root(std::string rootIN){root = rootIN;}

void MonteCarloHelper::print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt,
                                       std::string sys_name){
    
    unsigned long long sys_steps_global = 0;
    unsigned long long sys_steps_accpt_global = 0;
    
    MPI_Reduce(&sys_steps,&sys_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&sys_steps_accpt,&sys_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        std::cout << "\t System " + sys_name +" Acceptance Ratio: " << 100*(double)sys_steps_accpt_global/sys_steps_global << std::endl;
    }
}
void MonteCarloHelper::print_avg_energy(double estimator_total, unsigned long long num_steps){
    
    double estimator_total_global = 0;
    unsigned long long num_steps_global = 0;
    double standev = 0;
    double standev_global = 0;
    
    MPI_Reduce(&estimator_total,&estimator_total_global,1,MPI_DOUBLE,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&num_steps,&num_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        estimator_total_global = estimator_total_global/num_steps_global;
    }
    
    MPI_Bcast(&estimator_total_global,1,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    standev = pow(estimator_total_global - (estimator_total/num_steps),2);
    MPI_Reduce(&standev,&standev_global,1,MPI_DOUBLE,MPI_SUM,root_proc, MPI_COMM_WORLD);

    if (my_id == root_proc) {
        if (num_procs > 1) {
            standev_global = sqrt(standev_global/(num_procs-1));
            
            std::cout << "\t Average Energy: " << estimator_total_global << std::endl;
            std::cout << "\t Standard Deviation: " << standev_global << std::endl;
            std::cout << "\t Using " << num_procs << " processors." << std::endl;
        }
        else{
            std::cout << "\t Average Energy: " << estimator_total_global << std::endl;
        }
    }
}
void MonteCarloHelper::write_estimator(vector <double> estimator,int interval){
    
    int size = estimator.size();
    
    std::string fileName = root + "estimator.txt";
    
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
void MonteCarloHelper::write_PSV(int nuc_beads,vector<double> Q, std::string system){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "PSV" + system + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    for(int bead=0; bead<nuc_beads; bead++){
        myFile << Q(bead) << std::endl;
    }
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully saved PSV to Results." << std::endl;
    }
    
    myFile.close();
}
void MonteCarloHelper::read_PSV(int nuc_beads, vector<double> &Q, std::string system){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "PSV" + system + quickConvert.str();
    
    std::ifstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    for(int bead=0; bead<nuc_beads; bead++){
        myFile >> Q(bead);
    }
    
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully read PSV file from Results." << std::endl;
    }
}
void MonteCarloHelper::write_MC_data(double estimator_total, unsigned long long num_steps){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "mcData" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile << estimator_total << std::endl;
    myFile << num_steps << std::endl;
    
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully saved MC data to Results." << std::endl;
    }
}
void MonteCarloHelper::read_MC_data(double &estimator_total, unsigned long long &num_steps){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "mcData" + quickConvert.str();
    
    std::ifstream myFile;
    myFile.open(fileName.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile >> estimator_total;
    myFile >> num_steps;
    
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully read MC datat from Results." << std::endl;
    }
}

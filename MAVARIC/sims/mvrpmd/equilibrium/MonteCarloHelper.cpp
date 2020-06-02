#include "MonteCarloHelper.h"

MonteCarloHelper::MonteCarloHelper(std::string root, int my_id, int num_procs, int root_proc)
    :root(root), my_id(my_id),
     num_procs(num_procs), root_proc(root_proc)
{}
void MonteCarloHelper::print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt){
    
    unsigned long long sys_steps_global = 0;
    unsigned long long sys_steps_accpt_global = 0;
    
    MPI_Reduce(&sys_steps,&sys_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&sys_steps_accpt,&sys_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        std::cout << "\t System Acceptance Ratio: " << 100*(double)sys_steps_accpt_global/sys_steps_global << std::endl;
    }
}
void MonteCarloHelper::print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt){
   
    unsigned long long elec_steps_global = 0;
    unsigned long long elec_steps_accpt_global = 0;
    
    MPI_Reduce(&elec_steps,&elec_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&elec_steps_accpt,&elec_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        std::cout << "\t Electronic Acceptance Ratio: " << 100*(double)elec_steps_accpt_global/elec_steps_global << std::endl;
    }
}
void MonteCarloHelper::print_avg_energy(double estimator_total, double sgn_total){
    
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

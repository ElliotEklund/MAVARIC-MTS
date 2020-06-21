#include "SamplingHelper.h"

SamplingHelper::SamplingHelper(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs), root_proc(root_proc)
{}

void SamplingHelper::set_root(std::string rootIN){root = rootIN;}

void SamplingHelper::print_sys_accpt(unsigned long long sys_steps,
                                     unsigned long long sys_steps_accpt,
                                     std::string sys_name){

    unsigned long long sys_steps_global = 0;
    unsigned long long sys_steps_accpt_global = 0;

    MPI_Reduce(&sys_steps,&sys_steps_global,1,MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    MPI_Reduce(&sys_steps_accpt,&sys_steps_accpt_global,1,
               MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);

    if (my_id == root_proc) {
        std::cout << "\t System " + sys_name + " Acceptance Ratio: "
        << 100*(double)sys_steps_accpt_global/sys_steps_global << std::endl;
    }
}
void SamplingHelper::read_PSV(int nuc_beads,int elec_beads,int num_states,
                              vector<double> &Q, matrix<double> &x,matrix<double> &p){
    
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
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            myFile >> x(bead,state);
        }
    }
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            myFile >> p(bead,state);
        }
    }
    
    myFile.close();
    
    if (my_id == root_proc) {
        std::cout << "\t Successfully read PSV file." << std::endl;
    }
}

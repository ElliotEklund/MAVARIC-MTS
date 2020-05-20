#include "SamplingHelper.h"

SamplingHelper::SamplingHelper(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs), root_proc(root_proc)
{}

void SamplingHelper::set_root(std::string rootIN){root = rootIN;}

void SamplingHelper::print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt,
                                       std::string sys_name){
    
    unsigned long long sys_steps_global = 0;
    unsigned long long sys_steps_accpt_global = 0;
    
    MPI_Reduce(&sys_steps,&sys_steps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    MPI_Reduce(&sys_steps_accpt,&sys_steps_accpt_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,root_proc, MPI_COMM_WORLD);
    
    if (my_id == root_proc) {
        std::cout << "\t System " + sys_name +" Acceptance Ratio: " << 100*(double)sys_steps_accpt_global/sys_steps_global << std::endl;
    }
}

void SamplingHelper::read_PSV(int nuc_beads, vector<double> &Q, std::string system){
    
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

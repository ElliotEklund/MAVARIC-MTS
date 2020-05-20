#include "mpi_wrapper.hpp"

mpi_wrapper::mpi_wrapper(int num_procs, int my_id, int root_proc)
    :num_procs(num_procs),
     my_id(my_id),
     root_proc(root_proc)
{}
void mpi_wrapper::write_vector(const vector<double> &v,std::string fileName){
    
    int v_size = v.size(); //length of v
    int all_sizes [num_procs]; //collection of v_size across all procs
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(&v_size,1,MPI_INT,&all_sizes,1,
               MPI_INT,root_proc,MPI_COMM_WORLD);

    int total_size = 0;//sum of all v_size's

    if (my_id==root_proc) {
        for (int i=0; i<num_procs; i++) {
            total_size += all_sizes[i];
        }
    }

    double v_final[total_size]; //holds combined elements across all procs
    int displs[num_procs]; //buffer space for a given proc
    displs[0] = 0;

    for (int i=1; i<num_procs; i++) {
        displs[i] = displs[i-1] + all_sizes[i-1];
    }

    MPI_Gatherv(&v(0),v_size, MPI_DOUBLE,
                &v_final[0], all_sizes, displs, MPI_DOUBLE,
                root_proc, MPI_COMM_WORLD);

    if (my_id==root_proc) {
        std::ofstream myFile;
        myFile.open(fileName);
        if (!myFile.is_open()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
        }
        for (int i=0; i<total_size; i++) {
            myFile << v_final[i] << std::endl;
        }
    }
}

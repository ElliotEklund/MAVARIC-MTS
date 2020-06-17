#ifndef trajs_io_hpp
#define trajs_io_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <string>
#include <fstream>
#include <stdlib.h>
#include "mpi.h"

using namespace boost::numeric::ublas;

/* Read the first v_size lines of file and return them in a vector.
   Comment lines are skipped over in the file.  */

inline vector<double> load_file(std::string file, unsigned long long v_size){
    
    std::ifstream myStream;
    myStream.open(file.c_str());
    
    if(!myStream.is_open()) {
        std::cout << "ERROR: Could not open: " << file << std::endl;
    }
    
    /* Skip past file comment lines*/
    std::string line; //current line
    int num_lines = 0; //number of comment lines
    bool keep_going = true;
    
    while (keep_going) {
        std::getline(myStream,line);
        std::size_t found = line.find("#",0);
        if (found != std::string::npos) {
            //line is a comment line; it can be skipped
            num_lines +=1;
        }
        else{keep_going = false;}
    }
    
    myStream.close();
    
    myStream.open(file.c_str());
    std::string junk;
    
    /* Skip over comment lines */
    for (int i=0; i<num_lines; i++) {
        myStream >> junk;
    }
    
    /* Read in data*/
    vector<double> v(v_size,0);
    for (int i=0; i<v_size; i++){
        myStream >> v(i);
    }
    
    myStream.close();
    
    return v;
}

/* Distribute v_global over all processors. processor 0 will recieve the
first vl_size entries of v_global, processor 1 will recieve the next v1_size
elements, and so on. Note: vl_size maybe different for each processor,
The resulting local vector is returned. */
inline vector<double> dist_vec(vector<double> v_global,int vl_size,
                               int my_id, int num_procs,int root_proc){
    
    int all_sizes [num_procs]; //collection of v_size across all procs
    
    MPI_Gather(&vl_size,1,MPI_INT,&all_sizes,1,
               MPI_INT,root_proc,MPI_COMM_WORLD);
    
    int total_size = 0;//sum of all v_size's

    if (my_id==root_proc) {
        for (int i=0; i<num_procs; i++) {
            total_size += all_sizes[i];
        }
    }
    
    int displs[num_procs]; //buffer space for a given proc
    displs[0] = 0;

    for (int i=1; i<num_procs; i++) {
        displs[i] = displs[i-1] + all_sizes[i-1];
    }
    
    vector<double> v_local(vl_size,0.0);

    MPI_Scatterv(&v_global(0),all_sizes, displs, MPI_DOUBLE, &v_local(0),
                 vl_size, MPI_DOUBLE, root_proc, MPI_COMM_WORLD);
    
    return v_local;
}

/* Read the first vg_size lines from file and distribute them over all
 processors. Each process recieves a vector with those distributed values. */
inline vector<double> get_trajs(std::string file, unsigned long long vg_size,
                         int vl_size, int my_id, int num_procs,int root_proc){
    
    vector<double> v_global;
    if (my_id==root_proc) {
        v_global.resize(vg_size,0);
    }
    
    v_global = load_file(file,vg_size);
    vector<double> v_local(vl_size,0);
    v_local = dist_vec(v_global,vl_size,my_id,num_procs,root_proc);
    
    return v_local;
}

/* Return a reformated version of v. The reformatted vector will be a vector
 of length num_trajs. Each element of the vector is also a vector of length
 num_beads*/
inline vector<vector<double> > reformat(vector<double> &v,
                                        unsigned long long num_trajs,
                                        int num_beads){
    
    vector<vector<double> > a (num_trajs,zero_vector<double>(num_beads));
    int s = 0; //stride
    
    for (int traj=0; traj<num_trajs; traj++) {
        s = traj*num_beads;
        for (int bead=0; bead<num_beads; bead++) {
            a(traj)(bead) = v(s + bead);
        }
    }
    
    return a;
}

/* Return a reformated version of v. The reformatted vector will be a vector
of length num_trajs. Each element of the vector is a matrix of dimension
num_beads x num_states*/
inline vector<matrix<double> > reformat(vector<double> &v,
                                        unsigned long long num_trajs,
                                        int num_beads, int num_states){
    
    vector<matrix<double> > a (num_trajs,zero_matrix<double>(num_beads,num_states));
    int s1 = 0; //stride
    int s2 = 0; //stride
    
    for(int traj=0; traj<num_trajs; traj++){
        s1 = traj*num_beads*num_states;
        for(int bead=0; bead<num_beads; bead++){
            s2 = bead*num_states;
            for(int state=0; state<num_states; state++){
                a(traj)(bead,state) = v(s1 + s2 + state);
            }
        }
    }
    
    return a;
}

inline vector<vector<double> > get_trajs_reformat(std::string file,
                                                  unsigned long long vg_size,
                                                  int vl_size, int my_id,
                                                  int num_procs,
                                                  int root_proc, int num_trajs,
                                                  int num_beads){
    
    
    vector<double> v_trajs(vl_size,0);
    v_trajs = get_trajs(file, vg_size, vl_size,my_id,num_procs,root_proc);
    return reformat(v_trajs,num_trajs,num_beads);
}
inline vector<matrix<double> > get_trajs_reformat(std::string file,
                                                  unsigned long long vg_size,
                                                  int vl_size, int my_id,
                                                  int num_procs,
                                                  int root_proc, int num_trajs,
                                                  int num_beads,int num_states){
    
    vector<double> v_trajs(vl_size,0);
    v_trajs = get_trajs(file, vg_size, vl_size,my_id,num_procs,root_proc);
    return reformat(v_trajs,num_trajs,num_beads,num_states);
}

#endif

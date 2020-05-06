#include "aggregate.hpp"

aggregate::aggregate(){
    
}

void aggregate::add_calc(std::string name, int num_cols, int num_rows){
    
    matrix<double> * v = new matrix<double> (num_rows,num_cols + 1,0.0);
    myMap.insert(std::pair<std::string, matrix<double> * >(name,v));
}

void aggregate::collect(std::string name, int row, const vector<double> & v0,
                        const vector<double> & v, const double & sgnTheta){
    
    int num_cols = v.size();

    for (int col=0; col<num_cols; col++) {
        (*myMap[name])(row,col) += sgnTheta * v0(col) * v(col);
    }
    
    (*myMap[name])(row,num_cols) += sgnTheta;
}

void aggregate::merge_collections(int root_process,int my_id, std::string root){
    
    std::map<std::string, matrix<double> *>::iterator itr;
    
    std::cout << "merging" << std::endl;
    
    for (itr = myMap.begin(); itr != myMap.end(); itr++) {
        std::string name = itr->first;
        int num_rows = itr->second->size1();
        int num_cols = itr->second->size2();
        
        matrix<double> v_sum (num_rows,num_cols,0.0);
        
        MPI_Reduce(&(*itr->second)(0,0),&v_sum(0,0),num_rows*num_cols,
                        MPI_DOUBLE,MPI_SUM,root_process,MPI_COMM_WORLD);
        
        if (my_id==root_process) {
            std::string fileName = root + name;
            std::cout << fileName << std::endl;
            std::ofstream myFile;
            myFile.open(fileName);
            
            if (!myFile.is_open()) {
                std::cout << "ERROR: Could not open " << fileName << std::endl;
            }
            
            for (int row=0; row<num_rows; row++) {
                for (int col=0; col<num_cols; col++) {
                    myFile << v_sum(row,col) << " ";
                }
                myFile << std::endl;
            }
            myFile.close();
        }
    }
}

void aggregate::print_collection(std::string name){
    std::cout <<  (*myMap[name]) << std::endl;
}

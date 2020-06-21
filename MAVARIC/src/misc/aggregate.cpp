#include "aggregate.hpp"

aggregate::aggregate(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{}
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

void aggregate::collect(std::string name, int row, const double & v0,
                        const double & v, const double & sgnTheta){
    
    int num_cols = 1;
    
    for (int col=0; col<num_cols; col++) {
        (*myMap[name])(row,col) += sgnTheta * v0 * v;
    }
    
    (*myMap[name])(row,num_cols) += sgnTheta;
}

/* HACK ALERT!!!  I am currently hard ss (stride). This should
 be fixed at a later time
 num_trajs is the global num_trajs*/
void aggregate::merge_collections(int root_process,int my_id, std::string root,
                                  double dt, double ss, unsigned long long num_trajs){
        
    std::map<std::string, matrix<double> *>::iterator itr;
        
    for (itr = myMap.begin(); itr != myMap.end(); itr++) {
        std::string name = itr->first;
        int num_rows = itr->second->size1();
        int num_cols = itr->second->size2();

        double v[num_rows*num_cols];
        double v_sum[num_rows*num_cols];
        
        for (int col=0; col<num_cols; col++) {
            int stride = col*num_rows;
            for (int row=0; row<num_rows; row++){
                v[stride + row] = (*itr->second)(row,col);
            }
        }
        
        MPI_Reduce(&v[0],&v_sum[0],num_rows*num_cols,
                        MPI_DOUBLE,MPI_SUM,root_process,MPI_COMM_WORLD);
        
        if (my_id==root_process) {
            std::string fileName = root + name;
            std::ofstream myFile;
            myFile.open(fileName.c_str());
            
            if (!myFile.is_open()) {
                std::cout << "ERROR: Could not open " << fileName << std::endl;
            }
            
            myFile << "#dt:" << dt << std::endl;
            myFile << "#num_trajs:" << num_trajs << std::endl;

            
            int stride1 = 0;
            int stride2 = 0;
            stride2 = (num_cols-1)*num_rows;
            
            for (int row=0; row<num_rows; row++) {
                myFile << row*dt*ss << " ";
                for (int col=0; col<num_cols-1; col++) {
                    stride1 = col*num_rows;
                    myFile << v_sum[stride1 + row]/v_sum[stride2 + row] << " ";
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
void aggregate::write_errors(std::string name,int num_samples,int num_errors,
                             double dt,int stride,std::string output_dir){
    
    int num_cols = myMap[name]->size2();
    vector<vector<double> > errors (num_cols,zero_vector<double>(num_errors));
    vector<vector<double> > avgs (num_cols-1,zero_vector<double>(num_errors));
    vector<vector<double> > stdevs (num_cols-1,zero_vector<double>(num_errors));
    vector<vector<double> > temp (num_samples,zero_vector<double>(num_errors));
    vector<vector<double> > sgnTheta (num_samples,zero_vector<double>(num_errors));

    for (int i=0; i<num_cols; i++) {
        errors(i) = get_samples(name,i,num_errors);
    }

    sgnTheta = sum_data(errors(num_cols-1),num_samples,num_errors);
    
    for (int col=0; col<num_cols-1; col++) {
        temp = sum_data(errors(col),num_samples,num_errors);
        avgs(col) = average(temp,num_samples,num_errors,sgnTheta);
    }
    
    for (int col=0; col<num_cols-1; col++) {
        temp = sum_data(errors(col),num_samples,num_errors);
        stdevs(col) = stdev(temp,avgs(col),sgnTheta,num_samples,num_errors);
    }

    if (my_id == root_proc) {
        std::ofstream myStream;
        std::string file_name = output_dir + name + "_error";
        myStream.open(file_name.c_str());
        if (!myStream.is_open()) {
            std::cout << "ERROR: Could not open " << file_name << std::endl;
        }
        
        myStream << "#num_samples:" << num_samples << std::endl;
        
        int num_rows = myMap[name]->size1();
        int stride2 = floor(double(num_rows)/num_errors);
        
        for (int i=0; i<num_errors; i++) {
            myStream << i*dt*stride*stride2<< " ";
            for (int j=0; j<num_cols-1; j++) {
                myStream << stdevs(j)(i) << " ";
            }
            myStream << std::endl;
        }
        myStream.close();
    }
}
vector<double> aggregate::get_samples(std::string name,int col,int num_errors){
    
    int num_rows = myMap[name]->size1();
    int num_cols = myMap[name]->size2();

    int stride = floor(double(num_rows)/num_errors);
    vector<double> v_error(num_errors,0);
    
    int count = 0;
    for (int i=0; i<num_rows; i++) {
        if (i % stride ==0) {
            v_error(count) = (*myMap[name])(i,col);
            count += 1;
        }
    }
    return v_error;
}
vector<vector<double> > aggregate::sum_data(vector<double> v_error,
                                             int num_samples, int num_errors){
        
    if (num_samples > num_procs) {
        if (my_id==root_proc) {
            std::cout << "ERROR: num_procs must be greater than or equal to "
            "num_samples" << std::endl;
        }
    }
    
    /* Number of groups that can be made with num_procs each containining numsamples*/
    int num_groups = floor(double(num_procs)/num_samples);
    int r = num_procs % num_samples; //remainder procs

    vector<double> v_error_global (num_errors*num_procs);

    MPI_Allgather(&v_error(0),num_errors,MPI_DOUBLE,
                  &v_error_global(0),num_errors,MPI_DOUBLE,
                  MPI_COMM_WORLD);
    
    vector<vector<double> > v_error_sum(num_samples,zero_vector<double>(num_errors));
    
    for (int i=0; i<num_samples; i++) {
        for (int j=0; j<num_groups; j++) {
            for (int k=0; k<num_errors; k++) {
                v_error_sum(i)(k) += v_error_global(i*num_groups*num_errors + j*num_errors + k);
            }
        }
    }
    
    return v_error_sum;
}
vector<double> aggregate::average(vector<vector<double> > v, int num_samples,
                                  int num_errors,vector<vector<double> > denom){
    
    vector<double> avg_v (num_errors,0);
    vector<double> avg_denom (num_errors,0);

    for (int i=0; i<num_samples; i++) {
        for (int j=0; j<num_errors; j++) {
            avg_v(j) += v(i)(j);
            avg_denom(j) += denom(i)(j);
        }
    }
    
    for (int i=0; i<num_errors; i++) {
        avg_v(i) = avg_v(i)/avg_denom(i);
    }
    
    return avg_v;
}
vector<double> aggregate::stdev(vector<vector<double> > v,vector<double> avgs,
                                vector<vector<double> > denom, int num_samples,
                                int num_errors){
    
    vector<double> sigma (num_errors,0);
    
    for (int i=0; i<num_errors; i++) {
        for (int j=0; j<num_samples; j++) {
            double dif = avgs(i) - (v(j)(i)/denom(j)(i));
            sigma(i) += dif*dif;
        }
        sigma(i) = sqrt(sigma(i)/num_samples);
    }
    return sigma;
}

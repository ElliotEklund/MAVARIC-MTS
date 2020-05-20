#ifndef STEPPER_RPMD_B_H
#define STEPPER_RPMD_B_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "system_move.h"
#include "bath_move.h"

using namespace boost::numeric::ublas;

class Stepper_RPMD_B{
    
public:
    
private:
    
    system_move s_move;
    bath_move b_move;
    
    
};

#endif

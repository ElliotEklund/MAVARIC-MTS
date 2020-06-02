cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

set(ROOT_MAC "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC")

set(BOOST_MAC "/Users/ellioteklund/Downloads")

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

if(${MACHINE} STREQUAL "MAC")
    set(BOOST_DIR "${BOOST_MAC}/boost_1_72_0")
    set(HAMILTONIAN_DIR "${ROOT_MAC}/src/hamiltonian")
    set(MISC_DIR "${ROOT_MAC}/src/misc")
    set(CMAKE_CXX_FLAGS "-DMAC")
endif()

add_definitions( "-DNDEBUG -O3" )

# add the executable

if(${MACHINE} STREQUAL "MAC")
    add_executable(test test.cpp

        #hamiltonian directory
        ${HAMILTONIAN_DIR}/C_Matrix.cpp
        ${HAMILTONIAN_DIR}/M_Matrix_MTS.cpp
        ${HAMILTONIAN_DIR}/M_Matrix.cpp
        ${HAMILTONIAN_DIR}/StateDepPots.cpp
        ${HAMILTONIAN_DIR}/Theta_MTS.cpp
        ${HAMILTONIAN_DIR}/Theta.cpp

        #misc directory
        ${MISC_DIR}/MainHlpr.cpp
        ${MISC_DIR}/pop_estimators.cpp
        ${MISC_DIR}/aggregate.cpp
        ${MISC_DIR}/PopulationEstimator.cpp)
endif()

include_directories(${BOOST_DIR})
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)


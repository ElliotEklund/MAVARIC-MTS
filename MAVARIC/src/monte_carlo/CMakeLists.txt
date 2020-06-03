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
    set(MONTECARLO_DIR "${ROOT_MAC}/src/monte_carlo")
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
        ${HAMILTONIAN_DIR}/StateIndepPot.cpp
        ${HAMILTONIAN_DIR}/SpringEnergy.cpp
        ${HAMILTONIAN_DIR}/GTerm.cpp
        ${HAMILTONIAN_DIR}/Theta_MTS.cpp
        ${HAMILTONIAN_DIR}/Theta.cpp
        ${HAMILTONIAN_DIR}/theta_mixed.cpp
        ${HAMILTONIAN_DIR}/MVRPMD_MTS_Hamiltonian.cpp
        ${HAMILTONIAN_DIR}/mvrpmd_ham_mixed.cpp

	#Monte Carlo directory
        ${MONTECARLO_DIR}/system_step.cpp

        #misc directory
        ${MISC_DIR}/MainHlpr.cpp
        ${MISC_DIR}/pop_estimators.cpp
        ${MISC_DIR}/aggregate.cpp
        ${MISC_DIR}/PopulationEstimator.cpp)
endif()

include_directories(${BOOST_DIR})
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)
include_directories(${MONTECARLO_DIR}/Headers)


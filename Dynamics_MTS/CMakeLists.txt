cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

set(ROOT_MAC "/Users/ellioteklund/Desktop")
set(ROOT_LIN "/home/elliot/Desktop")

set(BOOST_MAC "/Users/ellioteklund/Downloads")
set(BOOST_LIN "/home/elliot")

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

if(${MACHINE} STREQUAL "MAC")
    set(INCLUDE_DIRS "${BOOST_MAC}/boost_1_72_0")
    set(MY_HEADERS "${ROOT_MAC}/Dynamics_MTS_git/Dynamics_MTS/DynamicsParts/Headers")
    set(MY_SOURCE "${ROOT_MAC}/Dynamics_MTS_git/Dynamics_MTS/DynamicsParts/Sources")

elseif(${MACHINE} STREQUAL "LIN")
    set(INCLUDE_DIRS "${BOOST_LIN}/boost_1_72_0")
    set(MY_HEADERS "${ROOT_LIN}/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/DynamicsParts/Headers")
    set(MY_SOURCE "${ROOT_LIN}/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/DynamicsParts/Sources")
endif()

add_definitions( "-DNDEBUG -O3" )

# add the executable
add_executable(main main.cpp
    ${MY_SOURCE}/dSpring_dQ.cpp
    ${MY_SOURCE}/dCdelec.cpp
    ${MY_SOURCE}/GTerm.cpp
    ${MY_SOURCE}/C_Matrix.cpp
    ${MY_SOURCE}/M_Matrix.cpp
    ${MY_SOURCE}/M_Matrix_MTS.cpp
    ${MY_SOURCE}/dM_Matrix_MTS_dQ.cpp
    ${MY_SOURCE}/dM_Matrix_dQ.cpp
    ${MY_SOURCE}/StateDepPots.cpp
    ${MY_SOURCE}/dStateDep_dQ.cpp
    ${MY_SOURCE}/dStateIndep_dQ.cpp
    ${MY_SOURCE}/StateIndepPot.cpp
    ${MY_SOURCE}/SpringEnergy.cpp
    ${MY_SOURCE}/Theta_MTS.cpp
    ${MY_SOURCE}/dTheta_MTS_dElec.cpp
    ${MY_SOURCE}/dTheta_MTS_dQ.cpp
    ${MY_SOURCE}/Forces_MTS.cpp
    ${MY_SOURCE}/MVRPMD_MTS_Hamiltonian.cpp
    ABM_MVRPMD.cpp
    RK4_MVRPMD.cpp
    Dynamics.cpp)
include_directories(${INCLUDE_DIRS})
include_directories(${MY_HEADERS})


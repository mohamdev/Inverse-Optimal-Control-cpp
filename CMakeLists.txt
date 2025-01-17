cmake_minimum_required(VERSION 3.15)
project(optim-cpp-client VERSION 0.1.0)

add_executable(optim-cpp "")

target_sources(optim-cpp
    PUBLIC
        apps/main_IOC_bayesopt.cpp
        include/mechanical_model.h
        include/utils.h
        include/cost-functions.h
        include/NLP_DOC_constr_humantable.h
        include/polynomial_splines.h
        include/bayesopt_ioc.h
        src/mechanical_model.cpp
        src/cost-functions.cpp
        src/NLP_DOC_constr_humantable.cpp
        src/polynomial_splines.cpp
        src/utils.cpp
)

find_library(PIN_LIB pinocchio) #pinocchio
find_library(HPP_FCL hpp-fcl) #pinocchio
find_library(URDFDOM_LIB urdfdom_model) #urdfdom
# find_library(BOOST_FS_LIB boost_filesystem) #boost
# find_library(BOOST_SER_LIB boost_serialization) #boost
# find_library(BOOST_SYS_LIB boost_system)
find_library(DL_LIB dl)
find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)
find_library(M_LIB m)
find_library(CPPAD_LIB cppad_lib)
find_library(COINMUMPS_LIB coinmumps)
find_library(OMNI_ORB omniORB4)
find_library(OMNI_THREAD omnithread3)
find_library(IPOPT_LIB ipopt)
# find_library(NDCURVES_LIB ndcurves)

include_directories(${OpenCV_INCLUDE_DIRS})

# This call requires CMake 3.13 or later, see next section
target_link_libraries(optim-cpp PUBLIC ${PIN_LIB} 
                    ${URDFDOM_LIB} 
                    ${DL_LIB} 
                    ${BLAS_LIB}
                    ${LAPACK_LIB} 
                    ${M_LIB}
                    ${CPPAD_LIB}
                    ${COINMUMPS_LIB}
                    ${OMNI_ORB}
                    ${IPOPT_LIB}
                    )


target_compile_options(optim-cpp PUBLIC -c 
                    -fmessage-length=0 
                    -DPINOCCHIO_WITH_URDFDOM=ON 
                    -DBOOST_MPL_LIMIT_LIST_SIZE=30
                    -DHAVE_CSTDDEF
                    -I/usr/local/include/coin-or
                    -I/usr/include/coin
                    -I/usr/local/include/coin-or/mumps
                    -I/usr/include/x86_64-linux-gnu
                    -I/usr/local//usr/local/include
                    -I/usr/include/x86_64-linux-gnu/python2.7 
                    -I/usr/include/python2.7
                    -I/opt/openrobots/include
                    )

target_link_libraries(optim-cpp PUBLIC /usr/lib/x86_64-linux-gnu/libomnithread.so.4)
target_link_libraries(optim-cpp PUBLIC /usr/local/lib/libbayesopt.a)
target_link_libraries(optim-cpp PUBLIC /usr/local/lib/libnlopt.a)
target_include_directories(optim-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(optim-cpp PUBLIC /usr/local/include/)
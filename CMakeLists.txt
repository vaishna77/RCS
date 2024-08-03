project (direct)
set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-13")
set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-13")
cmake_minimum_required (VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
set(GCC_NO_COMMON_BLOCKS "NO")

# Provide the input file that you want to solve for:
# This can also be passed as a command line argument:
if(NOT INPUT_FILE)
  #set (INPUT_FILE "test/CapacitanceExtraction.cpp")
  #set (INPUT_FILE "test/tuneFillInRank.cpp")
  #set (INPUT_FILE "test/svd_vs_aca.cpp")
  #set (INPUT_FILE "test/cube.cpp") # for RCS cube
  set (INPUT_FILE "test/sphere.cpp") # for RCS sphere
  #set (INPUT_FILE "test/debug_Matlab.cpp") # for 2D plate
  #set (INPUT_FILE "test/EFIE.cpp") # for 2D plate
  #set (INPUT_FILE "test/test.cpp")
endif()
# Here set the name for the output executable:
# This can also be passed as a command line argument:
if(NOT OUTPUT_EXECUTABLE)
    set (OUTPUT_EXECUTABLE "testEFIE")
endif()

# Set this to: "TWO", "THREE" #DIMENSION of the problem
# set (DIM "TWO")
 set (DIM "THREE")
# Set this to: "float", "double", "complex32" and "complex64"
#set (DTYPE "double")
set (DTYPE "complex64")

# Change this accordingly to point to the Eigen root folder
# By default, it will pick up the environment variable EIGEN_PATH
set(EIGEN_PATH "$ENV{EIGEN_PATH}")

################ FINDING MKL ####################

# Checking for MKL:
# Currently disabling MKL. Sometimes a linker error pops up
# TODO: Look into this
#if (DEFINED ENV{MKLROOT})
#    add_compile_definitions(MKL_ENABLED)
#    set (MKLROOT "$ENV{MKLROOT}")
#    set (LDFLAGS "${LDFLAGS} -DMKL_ENABLED -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group")
#    set (MKL_FLAGS "-DMKL_LP64 -DMKL_ENABLED -I${MKLROOT}/include")
#endif ()

################ FINDING EIGEN ####################
if(NOT EIGEN_PATH)
    message ( FATAL_ERROR "Please point the environment variable EIGEN_PATH to the root directory of your Eigen3 installation.")
endif()

################ SETTING DIMENSION TO USE:2 or 3 ####################
if(DIM STREQUAL "TWO")
    set(DIM_FLAG "-DUSE_TWO")
elseif(DIM STREQUAL "THREE")
    set(DIM_FLAG "-DUSE_THREE")
else()
    message(FATAL_ERROR "Invalid Choice for Dimension. Exiting...")
endif()

################ SETTING DATATYPE TO USE:FLOAT, DOUBLE AND COMPLEX ####################
if(DTYPE STREQUAL "float")
    set(DTYPE_FLAG "-DUSE_FLOAT")
elseif(DTYPE STREQUAL "double")
    set(DTYPE_FLAG "-DUSE_DOUBLE")
elseif(DTYPE STREQUAL "complex32")
    set(DTYPE_FLAG "-DUSE_COMPLEX32")
elseif(DTYPE STREQUAL "complex64")
    set(DTYPE_FLAG "-DUSE_COMPLEX64")
else()
    message(FATAL_ERROR "Invalid Choice for Datatype. Exiting...")
endif()

################ COMPILERS: APPROPRIATELY USE INTEL OR GNU ####################

# Toggle this to 1 when Intel Compilers are available:
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DIM_FLAG} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -O3 -ipo -qopenmp -xHost -funroll-loops -g -liomp5 -lpthread -ldl")
else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DIM_FLAG} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -g -fopenmp -lpthread -ldl")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DIM_FLAG} ${DTYPE_FLAG} ${MKL_FLAGS} --std=c++17 -O3 -funroll-loops -g -fopenmp -lpthread -ldl")
endif ()

# Flags for testing coverage information
# Need by coveralls:
set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS} --coverage")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS} --coverage")

################ ADDING OF THE NECESSARY HEADER FOLDERS ####################
# Adding folders which contain necessary header files:
include_directories("${EIGEN_PATH}"
                    "include"
                   )

################ MENTIONING C++ FILES TO LINK WITH HEADERS + FILE TO BUILD ####################

set(Platelib_SRCS
    ${PROJECT_SOURCE_DIR}/src/includes.cpp
    ${PROJECT_SOURCE_DIR}/src/kernel.cpp
    ${PROJECT_SOURCE_DIR}/src/Patch.cpp
    ${PROJECT_SOURCE_DIR}/src/EFIEMatrix.cpp
    #${PROJECT_SOURCE_DIR}/src/ElectrostaticsMatrix.cpp
    ${PROJECT_SOURCE_DIR}/src/AIFMMBox.cpp
    ${PROJECT_SOURCE_DIR}/src/ACA.cpp
    #${PROJECT_SOURCE_DIR}/src/basicKernel.cpp
    ${PROJECT_SOURCE_DIR}/src/AIFMM.cpp
    ${PROJECT_SOURCE_DIR}/src/rrqr.cpp
   )

if(NOT("${CMAKE_BUILD_TYPE}" STREQUAL "COVERAGE"))
    #add_subdirectory(examples)
    #add_subdirectory(test)
endif()

if(INPUT_FILE)
  add_executable(${OUTPUT_EXECUTABLE}
                 ${INPUT_FILE}
                 ${Platelib_SRCS}
                )
endif()

################ SUMMARY OF BUILD INFORMATION ################

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("Input File       : " ${INPUT_FILE})
message("Data Type        : " ${DTYPE})
message("CXX Compiler     : " ${CMAKE_CXX_COMPILER})
message("CXX_FLAGS        : " ${CMAKE_CXX_FLAGS})
message("Eigen dir        : " ${EIGEN_PATH})
# message("MKL dir          : " ${MKLROOT})
message("Build dir        : " ${CMAKE_BINARY_DIR})
message("PROJECT_SOURCE_DIR dir        : " ${PROJECT_SOURCE_DIR})
message("")

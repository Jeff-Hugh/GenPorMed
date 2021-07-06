# CMake 
cmake_minimum_required (VERSION 3.9)

# Project information
project (GeneratePorousMedia)
set (Demo_VERSION_MAJOR 0)
set (Demo_VERSION_MINOR 1)
set(EXECUTABLE_NAME "run")

# specify the C++ standard 
set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED True) 

INCLUDE_DIRECTORIES( ./)
aux_source_directory(./ DIR_SRCS)

FIND_PACKAGE(OpenMP)

set(Boost_INCLUDE_DIR /usr/share/boost_1_72_0)
set(Boost_LIBRARY_DIR /usr/share/boost_1_72_0/lib)
FIND_PACKAGE( Boost 1.60 COMPONENTS regex iostreams system REQUIRED )

if(Boost_FOUND)
    INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
    
    MESSAGE( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    MESSAGE( STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    MESSAGE( STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
    MESSAGE( STATUS "Boost_LIBRARY_DIR = ${Boost_LIBRARY_DIR}.")

    add_executable(${EXECUTABLE_NAME} ${DIR_SRCS})
    target_link_libraries(${EXECUTABLE_NAME} ${Boost_LIBRARIES})
    target_link_libraries(${EXECUTABLE_NAME} OpenMP::OpenMP_CXX)
endif()
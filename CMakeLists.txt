cmake_minimum_required(VERSION 2.8.12)
project(EM_C++)
include_directories(include/)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lrt -lpthread ")

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})





add_executable(GENIE_GxG src/FAME_dev.cpp src/genotype.cpp src/storage.cpp)

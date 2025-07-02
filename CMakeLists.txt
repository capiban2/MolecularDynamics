cmake_minimum_required(VERSION 3.14)
project(MPI_MD_Project LANGUAGES CXX)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(include)
include_directories(${MPI_INCLUDE_PATH})

add_executable(md_sim main.cpp)

target_link_libraries(md_sim ${MPI_CXX_LIBRARIES})

# Optionally, enable MPI compiler wrappers
if(MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

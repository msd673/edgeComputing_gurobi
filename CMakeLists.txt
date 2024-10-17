cmake_minimum_required(VERSION 3.5)
project(gurobi_EC)

link_directories(/data12/home/msd/gurobi/gurobi1001/linux64/lib)
include_directories(/data12/home/msd/gurobi/gurobi1001/linux64/include)

set(CMAKE_CXX_STANDARD 14)

add_executable(gurobi_EC gurobi_EC.cpp)
target_link_libraries(gurobi_EC libgurobi_g++5.2.a libgurobi100.so)

                            
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_FLAGS "-std=c++11 -g -lntl  -lgmp -fopenmp -lprimecount")

project(Weirds)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(HEADER_FILES opnalg.hpp)

add_executable(run
  main.cpp
  #subset_sum.cpp
	opanalg.cpp
  list.cpp
  tree.cpp
  tools.cpp
  expalg.cpp
)

target_link_libraries(run ntl m gmp primecount)

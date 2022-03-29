cmake_minimum_required (VERSION 3.10)
project(ZaVI_Bouldern VERSION 1.0.0
        DESCRIPTION "Boulder Repository"
        LANGUAGES CXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(json)
add_subdirectory(utils)
add_subdirectory(Particle_Filter)
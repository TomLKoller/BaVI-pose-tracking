cmake_minimum_required (VERSION 3.1)
project("Zavi-Particle_Filter" VERSION 1.0.0
                DESCRIPTION "Particle Filter to use with common state Descriptor"
                LANGUAGES CXX)

set(CMAKE_CXX_FLAGS_RELEASE "-O3  -march=native -Wall ")
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
add_definitions(-DNDEBUG)
endif()
include(GNUInstallDirs)


add_executable(${PROJECT_NAME}_Smooth MACOSX_BUNDLE realDataSmoothing.cpp)
#add_executable(${PROJECT_NAME}_Smooth MACOSX_BUNDLE checkOutlier.cpp)

  target_include_directories(${PROJECT_NAME}_Smooth PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> ../utils ../json/include)

target_compile_features(${PROJECT_NAME}_Smooth PUBLIC cxx_std_17)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Boost 1.37 REQUIRED COMPONENTS system)
find_package(ADEKF_VIZ CONFIG REQUIRED)
find_package(Ceres REQUIRED)
find_package(libconfig++ REQUIRED)


find_package(Qt5 COMPONENTS Core Gui Widgets PrintSupport Svg Xml  REQUIRED)

# include JKQTPlotter
find_package(JKQTCommonSharedLib REQUIRED)
find_package(JKQTMathTextSharedLib REQUIRED)
find_package(JKQTPlotterSharedLib REQUIRED)
#Include Threads to start ekf loop in other thread
find_package(Threads REQUIRED)
#find VTK
find_package(VTK REQUIRED)

# include all components


target_include_directories(${PROJECT_NAME}_Smooth PUBLIC ${EIGEN3_INCLUDE_DIR} ~/repositories/ADEKF/include/ADEKF/ ../BatchViterbi/include ~/repositories/ADEKF/ADEKF_VIZ/include)
target_link_libraries(${PROJECT_NAME}_Smooth PUBLIC utils ${Boost_LIBRARIES}  nlohmann_json::nlohmann_json ADEKF_VIZ::ADEKF_VIZ ${CMAKE_THREAD_LIBS_INIT} ceres  ${LIBCONFIG++_LIBRARIES})


target_compile_definitions(${PROJECT_NAME}_Smooth PUBLIC MANIFOLD_WITH_CERES)






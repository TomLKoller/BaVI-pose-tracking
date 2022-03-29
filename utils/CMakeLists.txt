cmake_minimum_required(VERSION 3.0)
project(utils)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

find_package(glog REQUIRED)
find_package(libconfig++ REQUIRED COMPONENTS Setting Config)
#SET(BOOST_DIR "C:/Users/Public/Arbeit/boost_1_55_0" )
#SET(EIGEN_DIR "C:/Users/Public/Arbeit/eigen-3.3.7")
# Specify critical OpenSceneGraph build variables.
# Reference OpenSceneGraph includes and libraries.
#INCLUDE_DIRECTORIES(${BOOST_DIR})
#INCLUDE_DIRECTORIES(${EIGEN_DIR})
#INCLUDE_DIRECTORIES(${OSG_BUILD_DIR}/include)

add_library(utils SHARED
        ZaVI_Utils.cpp
        ZaviConfig.cpp
        threaded_objects.cpp
        CSV_Reader.cpp
        datasetFormat.cpp
        ./Plugins/ARTPlugin.cpp
        ./Plugins/FileIMUPlugin.cpp
        ./Plugins/ResultLoader.cpp
        ./oriestimu_src/CsgOriEstIMU.cpp 
        ./oriestimu_src/CSGutilities.cpp
        )

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
add_definitions(-DNDEBUG)
endif()

set_target_properties(utils PROPERTIES PUBLIC_HEADER *.h)


FIND_PACKAGE(OpenSceneGraph)
find_package(ADEKF_VIZ CONFIG REQUIRED)
find_package(Qt5 COMPONENTS Core Gui Widgets PrintSupport Svg Xml  REQUIRED)

# include JKQTPlotter
find_package(JKQTCommonSharedLib REQUIRED)
find_package(JKQTMathTextSharedLib REQUIRED)
find_package(JKQTPlotterSharedLib REQUIRED)
#Include Threads to start ekf loop in other thread
find_package(Threads REQUIRED)
find_package(Eigen3)
#find VTK
find_package(VTK COMPONENTS
        vtkCommonColor
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersSources
        vtkInteractionStyle
        vtkRenderingCore
        vtkRenderingFreeType
        vtkRenderingContextOpenGL2 QUIET)
find_package(Boost REQUIRED COMPONENTS system)

target_include_directories(utils PUBLIC . ../json/include ${EIGEN3_INCLUDE_DIR} ~/repositories/ADEKF/include/ADEKF)
TARGET_LINK_LIBRARIES(utils  ${VTK_LIBRARIES}  ADEKF_VIZ::ADEKF_VIZ  nlohmann_json::nlohmann_json stdc++fs  ${Boost_LIBRARIES} libconfig::config++)
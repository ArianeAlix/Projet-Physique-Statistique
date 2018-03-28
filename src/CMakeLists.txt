project(modele)

cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(POLICY CMP0043)
    cmake_policy( SET CMP0043  OLD )
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Add specific Find.cmake modules 
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ) 

find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

find_package(Qt5 QUIET COMPONENTS Xml OpenGL)
find_package(OpenGL)


add_definitions(-DBOOST_NO_CXX11_TEMPLATE_ALIASES)


# The executable itself
add_executable( modele main.cpp)


SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Qt/5.10.1/msvc2017_64/lib")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Qt/5.10.1/msvc2017_64/bin")
SET(CMAKE_INCLUDE_DIRECTORIES_BEFORE "C:/Qt/5.10.1/msvc2017_64/include")

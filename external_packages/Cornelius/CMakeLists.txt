cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(libCornelius)

# for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set(SOURCES
    GeneralGeometryElement.cpp
    Line.cpp
    Polygon.cpp
    Polyhedron.cpp
    Square.cpp
    Cube.cpp
    Hypercube.cpp
    Cornelius.cpp)

add_library(Cornelius SHARED ${SOURCES})

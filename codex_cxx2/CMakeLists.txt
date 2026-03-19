cmake_minimum_required(VERSION 3.16)
project(codex_cxx2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CODEX_USE_MPI "Enable MPI parallel tracing over particles" OFF)

find_package(NetCDF REQUIRED)
if (CODEX_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

find_path(NETCDF_C_INCLUDE_DIR
    NAMES netcdf.h
    HINTS
        /opt/homebrew/include
        /usr/local/include)

find_library(NETCDF_C_LIBRARY
    NAMES netcdf
    HINTS
        /opt/homebrew/lib
        /usr/local/lib)

if (NOT NETCDF_C_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find NetCDF C headers (netcdf.h)")
endif()

if (NOT NETCDF_C_LIBRARY)
    message(FATAL_ERROR "Could not find NetCDF C library (libnetcdf)")
endif()

add_executable(trace_fieldlines
    src/main.cpp
    src/Interpolator.cpp
    src/AparData.cpp
    src/AparFieldModel.cpp
    src/FieldLineIntegrator.cpp
    src/PunctureDetector.cpp
    src/PoincareOutput.cpp
    src/MatlabComparator.cpp
    src/ValidationSuite.cpp)

target_include_directories(trace_fieldlines
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${NETCDF_C_INCLUDE_DIR}
        ${NetCDF_INCLUDE_DIRS})

target_link_libraries(trace_fieldlines
    PRIVATE
        ${NETCDF_C_LIBRARY}
        ${NetCDF_LIBRARIES})

if (CODEX_USE_MPI)
    target_compile_definitions(trace_fieldlines PRIVATE CODEX_USE_MPI=1)
    target_link_libraries(trace_fieldlines PRIVATE MPI::MPI_CXX)
endif()

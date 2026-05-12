cmake_minimum_required(VERSION 3.16)
project(codex_cxx2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CODEX_USE_MPI "Enable MPI parallel tracing over particles" OFF)
option(CODEX_USE_VISKORES "Enable Viskores runtime integration" OFF)

set(NETCDF_ROOT "" CACHE PATH "NetCDF C installation prefix")
if (CODEX_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()
if (CODEX_USE_VISKORES)
    find_package(Viskores REQUIRED CONFIG)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    find_package(Kokkos CONFIG QUIET)
    if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND TARGET viskores::openmp)
        # Some Viskores installs export -fopenmp=libomp, while this local
        # AppleClang build uses the flags detected by FindOpenMP instead.
        # Also avoid linking a second OpenMP runtime beside Viskores' runtime.
        set_target_properties(viskores::openmp PROPERTIES INTERFACE_COMPILE_OPTIONS "")
        set_target_properties(viskores::openmp PROPERTIES INTERFACE_LINK_LIBRARIES "")
    endif()
endif()

find_path(NETCDF_C_INCLUDE_DIR
    NAMES netcdf.h
    HINTS
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        /opt/homebrew/include
        /usr/local/include
    PATH_SUFFIXES include)

find_library(NETCDF_C_LIBRARY
    NAMES netcdf
    HINTS
        ${NETCDF_ROOT}
        $ENV{NETCDF_ROOT}
        /opt/homebrew/lib
        /usr/local/lib
    PATH_SUFFIXES lib lib64)

if (NOT NETCDF_C_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find NetCDF C headers (netcdf.h)")
endif()

if (NOT NETCDF_C_LIBRARY)
    message(FATAL_ERROR "Could not find NetCDF C library (libnetcdf)")
endif()

message(STATUS "NetCDF include dir: ${NETCDF_C_INCLUDE_DIR}")
message(STATUS "NetCDF library    : ${NETCDF_C_LIBRARY}")

add_executable(trace_fieldlines
    src/main.cpp
    src/Interpolator.cpp
    src/AparData.cpp
    src/AparFieldModel.cpp
    src/FieldLineIntegrator.cpp
    src/TracePostProcessor.cpp
    src/ViskoresFieldLineTracer.cpp
    src/PunctureDetector.cpp
    src/PoincareOutput.cpp
    src/MatlabComparator.cpp
    src/ValidationSuite.cpp)

target_include_directories(trace_fieldlines
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${NETCDF_C_INCLUDE_DIR}
)

target_link_libraries(trace_fieldlines
    PRIVATE
        ${NETCDF_C_LIBRARY}
)

if (CODEX_USE_MPI)
    target_compile_definitions(trace_fieldlines PRIVATE CODEX_USE_MPI=1)
    target_link_libraries(trace_fieldlines PRIVATE MPI::MPI_CXX)
endif()

if (CODEX_USE_VISKORES)
    target_compile_definitions(trace_fieldlines PRIVATE CODEX_USE_VISKORES=1)
    target_compile_options(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_OPTIONS>)
    target_compile_definitions(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_DEFINITIONS>)
    target_include_directories(trace_fieldlines PRIVATE $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_INCLUDE_DIRECTORIES>)
    # Avoid loading two different OpenMP runtimes in one process.
    # Viskores already links against the OpenMP runtime from its own build.
    target_link_libraries(trace_fieldlines PRIVATE viskores::cont viskores::worklet)
    if (TARGET Kokkos::kokkos)
        target_link_libraries(trace_fieldlines PRIVATE Kokkos::kokkos)
    endif()
    viskores_add_target_information(trace_fieldlines
        DROP_UNUSED_SYMBOLS
        MODIFY_CUDA_FLAGS
        DEVICE_SOURCES
            src/ViskoresFieldLineTracer.cpp)
endif()

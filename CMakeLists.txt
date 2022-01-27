# minimum version of cmake
cmake_minimum_required(VERSION 3.10)

# project name
project(EIC-Impact-Studies CXX)

# set cpp standard
set(CMAKE_CXX_STANDARD 17)

# add compilator options
add_compile_options()

# add linking options
if (UNIX)
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
		set(UNIX_LIBRARIES stdc++fs)
	endif()
endif ()

# additional modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# find ROOT
find_package(ROOT COMPONENTS Core Hist REQUIRED)

# find HepMC3
find_package(HepMC3 REQUIRED)

# includes
include_directories("${PROJECT_SOURCE_DIR}/include" ${ROOT_INCLUDE_DIRS} ${HEPMC3_INCLUDE_DIR})

# src
file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# executable
add_executable(eic_impact_studies "${SOURCES}")
set_target_properties(eic_impact_studies PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

target_link_libraries(eic_impact_studies "${ROOT_LIBRARIES}" "${HEPMC3_LIBRARIES}" "${UNIX_LIBRARIES}")

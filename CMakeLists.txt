cmake_minimum_required(VERSION 3.18)
project(radiatedPower VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate the ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Locate Boost libraries
find_package(Boost 1.73.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Add the FFTW ROOT library
find_library(ROOTFFTW_LIB RootFftwWrapper REQUIRED)
message("ROOTFFTW_LIB is set to ${ROOTFFTW_LIB}")

# Add libraries
add_subdirectory(BasicFunctions)
list(APPEND EXTRA_LIBS BasicFunctions)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/BasicFunctions")

add_subdirectory(FieldClasses)
list(APPEND EXTRA_LIBS FieldClasses)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/FieldClasses")

add_subdirectory(SignalProcessing)
list(APPEND EXTRA_LIBS SignalProcessing)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/SignalProcessing")

add_subdirectory(Antennas)
list(APPEND EXTRA_LIBS Antennas)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/Antennas")

add_subdirectory(SingleElectronRadiation)
list(APPEND EXTRA_LIBS SingleElectronRadiation)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/SingleElectronRadiation")

add_subdirectory(ElectronDynamics)

add_subdirectory(executables)


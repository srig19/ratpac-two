cmake_minimum_required(VERSION 3.11...3.14)
include(CMakePackageConfigHelpers)
# Automatically get version of git tags.
execute_process(
    COMMAND git --git-dir=${CMAKE_SOURCE_DIR}/.git describe --tags --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE RATVERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND git --git-dir=${CMAKE_SOURCE_DIR}/.git describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE RATVERSION_SEMVERS
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT RATVERSION_SEMVERS)
  set(RATVERSION_SEMVERS 3.2.0)
endif()

project(ratpac
  DESCRIPTION "RAT, Plus Additional Code, Version 2"
  HOMEPAGE_URL "https://ratpac.readthedocs.io/en/latest/"
  VERSION ${RATVERSION_SEMVERS}
  LANGUAGES C CXX
)

# Fallback version:
if(NOT RATVERSION)
  set(RATVERSION ${PROJECT_VERSION})
endif()

###########################################################
# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Use RPATH for install
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
###########################################################
# Set default build type
set(default_build_type "Release")
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS -pg)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
add_compile_options(-Wall)
# Not turning on errors yet to ease transition
# add_compile_options(-Werror)

###########################################################
# Set CMAKE standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_OSX_ARCHITECTURES "x86_64")

###########################################################
# External libraries and dependencies
set(CRY_Enabled 0)
if(DEFINED ENV{CRYLIB})
  message(STATUS "Compiling with CRY enabled")
  set(CRY_Enabled 1)
  set(CRYLIBDIR $ENV{CRYLIB})
  set(CRYINCLUDE $ENV{CRYINCLUDE})
  set(CRYDATA $ENV{CRYDATA})
  set(CRY_LIBRARIES CRY)
  include_directories(${CRYINCLUDE})
  link_directories(${CRYLIBDIR})
endif()

set(TENSORFLOW_Enabled 0)
find_library(TENSORFLOW_LIBRARIES tensorflow)
if(TENSORFLOW_LIBRARIES)
  message(STATUS "Compiling with Tensorflow")
  set(TENSORFLOW_Enabled 1)
else()
  message(WARNING "Tensorflow Not Found")
  set(TENSORFLOW_LIBRARIES "")
endif()


###########################################################
# Options
set(RATPAC_BINARY_DIR ${PROJECT_BINARY_DIR}/bin CACHE PATH "Where to output executables")
set(RATPAC_LIBRARY_DIR ${PROJECT_BINARY_DIR}/lib CACHE PATH "Where to output library files")
set(RATPAC_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include CACHE PATH "Where to collect headers")
set(RATPAC_CMAKE_DIR ${PROJECT_BINARY_DIR}/cmake CACHE PATH "Where to cmake configuration")
option(RATPAC_ROOT "build bootstrapped root executable" ON)
option(RATPAC_RAT "build rat executable" ON)

###########################################################
# Set up third party libraries

if(ROOT_DIR)
  find_package(ROOT NO_DEFAULT_PATH CONFIG REQUIRED COMPONENTS Minuit2 ROOTTPython)
else()
  find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit2 ROOTTPython)
endif()
include(${ROOT_USE_FILE})

find_package(Geant4 REQUIRED)
find_package(Geant4 QUIET OPTIONAL_COMPONENTS qt)
include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIRS})
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${XercesC_INCLUDE_DIRS})

find_package(Threads REQUIRED)


###########################################################
# Configure header with our version information
configure_file(${CMAKE_SOURCE_DIR}/config/Config.hh.in
  ${CMAKE_SOURCE_DIR}/src/core/include/RAT/Config.hh
  @ONLY NEWLINE_STYLE UNIX)

# Configure local scripts
set(RATROOT ${CMAKE_INSTALL_PREFIX})
set(RATSHARE ${PROJECT_SOURCE_DIR})
configure_file(${CMAKE_SOURCE_DIR}/config/ratpac.sh.in
  ${PROJECT_SOURCE_DIR}/ratpac.sh
  @ONLY NEWLINE_STYLE UNIX)
configure_file(${CMAKE_SOURCE_DIR}/config/ratpac.csh.in
  ${PROJECT_SOURCE_DIR}/ratpac.csh
  @ONLY NEWLINE_STYLE UNIX)
#configure_file(${CMAKE_SOURCE_DIR}/config/RatpacConfig.cmake.in
#  ${PROJECT_SOURCE_DIR}/RatpacConfig.cmake 
#  @ONLY NEWLINE_STYLE UNIX)

###########################################################
# Compile our libraries and executables

add_subdirectory(src)
add_subdirectory(tools)

###########################################################
# Installation

# Configure global scripts
set(RATSHARE ${CMAKE_INSTALL_PREFIX}/share/RAT)
configure_file(${CMAKE_SOURCE_DIR}/config/ratpac.sh.in
  ${RATPAC_BINARY_DIR}/ratpac.sh
  @ONLY NEWLINE_STYLE UNIX)
configure_file(${CMAKE_SOURCE_DIR}/config/ratpac.csh.in
  ${RATPAC_BINARY_DIR}/ratpac.csh
  @ONLY NEWLINE_STYLE UNIX)
configure_file(${CMAKE_SOURCE_DIR}/config/RatpacConfig.cmake.in
  ${RATPAC_CMAKE_DIR}/RatpacConfig.cmake 
  @ONLY NEWLINE_STYLE UNIX)
write_basic_package_version_file(
  "${RATPAC_CMAKE_DIR}/RatpacConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion
)

# Install the data files
install(DIRECTORY ratdb/ DESTINATION ${RATSHARE}/ratdb
  PATTERN "ratdb/*")

install(DIRECTORY models/ DESTINATION ${RATSHARE}/models
  PATTERN "models/*")

# Install the python files
install(DIRECTORY python/ DESTINATION ${RATSHARE}/python
  PATTERN "python/*")

# Install rattest executable
install(PROGRAMS python/rattest.py DESTINATION bin RENAME rattest)

# Install the macro files
install(DIRECTORY macros/ DESTINATION ${RATSHARE}/macros
  PATTERN "macros/*")

# Install the doc files
install(DIRECTORY doc/ DESTINATION ${RATSHARE}/doc
  PATTERN "doc/*")

# Install env files
install(FILES ${RATPAC_BINARY_DIR}/ratpac.sh
  ${RATPAC_BINARY_DIR}/ratpac.csh
  DESTINATION bin)

install(FILES ${RATPAC_CMAKE_DIR}/RatpacConfig.cmake
  DESTINATION cmake)

install(FILES ${RATPAC_CMAKE_DIR}/RatpacConfigVersion.cmake
  DESTINATION cmake)

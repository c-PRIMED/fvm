#
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
PROJECT(Xdmf)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
#    cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${Xdmf_SOURCE_DIR}/CMake)

OPTION(BUILD_SHARED_LIBS "Build Shared XDMF Library" ON)
OPTION(XDMF_WRAP_PYTHON "Build a python wrapping interface" ON)

# We need ansi c-flags, especially on HP
SET(CMAKE_C_FLAGS "${CMAKE_ANSI_CFLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_REQUIRED_FLAGS ${CMAKE_ANSI_CFLAGS})

# Platform configuration tests.
INCLUDE(${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

SET(XDMF_INSTALL_DIR_CONFIG ${Xdmf_BINARY_DIR})
INCLUDE_DIRECTORIES(
  ${Xdmf_SOURCE_DIR}/libsrc
  ${Xdmf_BINARY_DIR}/libsrc)

SET (LIBRARY_OUTPUT_PATH ${Xdmf_BINARY_DIR}/bin/ CACHE PATH 
  "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${Xdmf_BINARY_DIR}/bin/ CACHE PATH 
  "Single output directory for building all executables.")

INCLUDE (${CMAKE_ROOT}/Modules/CheckLibraryExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)

SET(VTK_INSTALL_LIB_DIR "/lib")

#############################################################################
# VTK READER/WRITER to use .xdmf files in vtk
# and
# ZLIB library
# and
# XML2 library

#used by vtklibxml2
INCLUDE(vtkThirdParty)

OPTION(XDMF_BUILD_VTK "Build VTK reader" OFF)
IF(XDMF_BUILD_VTK)

  FIND_PACKAGE(VTK REQUIRED)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ENDIF(VTK_FOUND)

  # ZLIB Library
  SET(XDMF_ZLIB_LIBRARIES ${VTK_ZLIB_LIBRARIES})
  SET(HDF5_ZLIB_INCLUDE_DIRS)
  IF(NOT VTK_USE_SYSTEM_ZLIB)
    SET(HDF5_ZLIB_HEADER "vtkzlib/zlib.h")
  ENDIF(NOT VTK_USE_SYSTEM_ZLIB)

ELSE(XDMF_BUILD_VTK)

  SET(VTK_DIR)

  # ZLIB Library
  OPTION(XDMF_USE_SYSTEM_ZLIB "Use system Zlib" OFF)
  IF(XDMF_USE_SYSTEM_ZLIB)
    FIND_PACKAGE(ZLIB)
    SET(XDMF_ZLIB_LIBRARIES ${ZLIB_LIBRARY})
    SET(HDF5_ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR})
  ELSE(XDMF_USE_SYSTEM_ZLIB)
    SUBDIRS(Utilities/vtkzlib)
    SET(HDF5_ZLIB_INCLUDE_DIRS 
      "${Xdmf_BINARY_DIR}/Utilities;${Xdmf_SOURCE_DIR}/Utilities")
    INCLUDE_DIRECTORIES(${HDF5_ZLIB_INCLUDE_DIRS})
    SET(XDMF_ZLIB_LIBRARIES vtkzlib)
    SET(HDF5_ZLIB_HEADER "vtkzlib/zlib.h")
  ENDIF(XDMF_USE_SYSTEM_ZLIB)

ENDIF(XDMF_BUILD_VTK)

#############################################################################
# XML2 LIBRARY
OPTION(XDMF_SYSTEM_LIBXML2 "Use system libxml2" OFF)
IF(NOT XDMF_SYSTEM_LIBXML2)

  # Build our own vtklibxml2.
  SUBDIRS(Utilities/vtklibxml2)
  SET(LIBXML2_LIBRARY vtklibxml2)
  INCLUDE_DIRECTORIES(${Xdmf_BINARY_DIR}/Utilities/vtklibxml2)
  SET(XDMF_LIBXML2_INCLUDE_CONFIG  ${Xdmf_BINARY_DIR}/Utilities/vtklibxml2)
  SET(XDMF_LIBXML2_LIBRARIES_CONFIG  ${LIBXML2_LIBRARY})

ELSE(NOT XDMF_SYSTEM_LIBXML2)

  # Use a libxml2 that is installed on the machine somewhere.
  FIND_LIBRARY(LIBXML2_LIBRARY
    xml2
    /usr/lib64
    /usr/lib
    /usr/local/lib)
  FIND_PATH(LIBXML2_INCLUDE_PATH
    libxml
    /usr/include/libxml2
    /usr/include)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_PATH}
    ${LIBXML2_INCLUDE_PATH}/libxml2)
  SET(LIBXML2_LIBRARY ${LIBXML2_LIBS})
  SET(XDMF_LIBXML2_INCLUDE_CONFIG  ${LIBXML2_INCLUDE_PATH};${LIBXML2_INCLUDE_PATH}/libxml2)
  SET(XDMF_LIBXML2_LIBRARIES_CONFIG  ${LIBXML2_LIBRARY})

ENDIF(NOT XDMF_SYSTEM_LIBXML2)
IF(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE
    -D_CRT_NONSTDC_NO_DEPRECATE 
    -D_CRT_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
ENDIF(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400) 

#############################################################################
# MPI LIBRARY
OPTION(XDMF_BUILD_MPI "Build MPI Objects" ON)
SET(XDMF_MPI_LIBRARIES "")
IF(XDMF_BUILD_MPI)
    INCLUDE (${CMAKE_ROOT}/Modules/FindMPI.cmake)
    IF (MPI_INCLUDE_PATH)
      SET(XDMF_MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    ENDIF(MPI_INCLUDE_PATH)
    IF (MPI_LIBRARY)
      SET(XDMF_MPI_LIBRARIES ${XDMF_MPI_LIBRARIES} ${MPI_LIBRARY})
    ELSE (MPI_LIBRARY)
      MESSAGE(FATAL_ERROR 
        "Could not find the required MPI libraries, install MPI, or set MPI_INCLUDE_PATH and MPI_LIBRARY and MPI_EXTRA_LIBRARY, or turn XDMF_BUILD_MPI to OFF")
    ENDIF (MPI_LIBRARY)
    IF (MPI_EXTRA_LIBRARY)
      SET(XDMF_MPI_LIBRARIES ${XDMF_MPI_LIBRARIES} ${MPI_EXTRA_LIBRARY})
    ENDIF (MPI_EXTRA_LIBRARY)
ELSE(XDMF_BUILD_MPI)
    ADD_DEFINITIONS("-DXDMF_NO_MPI")
ENDIF(XDMF_BUILD_MPI)

#############################################################################
# HDF5 LIBRARY
OPTION(XDMF_SYSTEM_HDF5 "Use system HDF5" OFF)
IF(XDMF_SYSTEM_HDF5)
  OPTION(XDMF_SYSTEM_HDF5_IS_PARALLEL "HDF5 Built for MPI" OFF)
  IF (XDMF_SYSTEM_HDF5_IS_PARALLEL)
  ENDIF (XDMF_SYSTEM_HDF5_IS_PARALLEL)
  FIND_LIBRARY(HDF5_LIBRARY
    hdf5
    /usr/lib
    /opt/lib
    /usr/local/lib)
  FIND_PATH(HDF5_INCLUDE_PATH
    H5Ipublic.h
    /usr/include
    /opt/include
    /usr/local/include)
  MESSAGE(STATUS "Using system HDF5")
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_PATH})
  SET(HDF5_LIBRARY ${HDF5_LIBRARY})
  SET(XDMF_HDF5_INCLUDE_CONFIG  ${HDF5_INCLUDE_PATH})
  SET(XDMF_HDF5_LIBRARY_CONFIG  ${HDF5_LIBRARY})
ELSE(XDMF_SYSTEM_HDF5)
  INCLUDE_DIRECTORIES(${Xdmf_SOURCE_DIR}/Utilities/hdf5
    ${Xdmf_BINARY_DIR}/Utilities/hdf5)
  SUBDIRS(Utilities)
  SET(HDF5_LIBRARY vtkhdf5)
  SET(HDF5_CONFIG "${Xdmf_BINARY_DIR}/Utilities/hdf5/HDF5Config.cmake")
  IF(EXISTS "${HDF5_CONFIG}")
    INCLUDE("${HDF5_CONFIG}")
    ADD_DEFINITIONS(${HDF5_EXTRA_FLAGS})
  ENDIF(EXISTS "${HDF5_CONFIG}")
  SET(XDMF_HDF5_INCLUDE_CONFIG ${Xdmf_BINARY_DIR}/Utilities/hdf5;${Xdmf_SOURCE_DIR}/Utilities/hdf5)
  SET(XDMF_HDF5_LIBRARY_CONFIG  ${HDF5_LIBRARY})
ENDIF(XDMF_SYSTEM_HDF5)

#############################################################################
# stdlib and std namespace
IF(XDMF_BUILD_VTK)
  # if using VTK, then copy the VTK value for VTK_USE_ANSI_STDLIB
  SET(XDMF_USE_ANSI_STDLIB "${VTK_USE_ANSI_STDLIB}" CACHE INTERNAL
    "ansi c++ stdlib flag")
ELSE(XDMF_BUILD_VTK)
  # if not using VTK, then do the test for ansi stream headers
  # here
  INCLUDE(TestForANSIStreamHeaders)
  IF(CMAKE_NO_ANSI_STREAM_HEADERS) 
    SET(XDMF_USE_ANSI_STDLIB FALSE CACHE INTERNAL "ansi c++ stdlib flag")
  ELSE(CMAKE_NO_ANSI_STREAM_HEADERS)
    SET(XDMF_USE_ANSI_STDLIB TRUE CACHE INTERNAL "ansi c++ stdlib flag")
  ENDIF(CMAKE_NO_ANSI_STREAM_HEADERS)
ENDIF(XDMF_BUILD_VTK)
INCLUDE(TestForSTDNamespace)

# Enforce 0/1 as only possible values.  Needed for the means by which
# Configure.hxx is implemented.
IF(CMAKE_NO_STD_NAMESPACE)
  SET(XDMF_NO_STD_NAMESPACE 1)
ELSE(CMAKE_NO_STD_NAMESPACE)
  SET(XDMF_NO_STD_NAMESPACE 0)
ENDIF(CMAKE_NO_STD_NAMESPACE)

#############################################################################
# Now build xdmf library
SUBDIRS(libsrc)

IF(XDMF_BUILD_VTK)
  SUBDIRS(vtk)
  INCLUDE_DIRECTORIES(${Xdmf_SOURCE_DIR}/vtk 
    ${Xdmf_BINARY_DIR}/vtk)
ENDIF(XDMF_BUILD_VTK)

#############################################################################
# Testing
ENABLE_TESTING()
CONFIGURE_FILE(${Xdmf_SOURCE_DIR}/CTestCustom.cmake.in 
  ${Xdmf_BINARY_DIR}/CTestCustom.cmake @ONLY)
INCLUDE(CTest)
OPTION(VTK_USE_DISPLAY 
  "Turn this option off and tests will not popup windows" ON)
MARK_AS_ADVANCED(VTK_USE_DISPLAY)

#############################################################################
# Examples
OPTION(XDMF_BUILD_EXAMPLES "Build C++ Examples" OFF)
IF(XDMF_BUILD_EXAMPLES)
  SUBDIRS(Examples/Cxx)
ENDIF(XDMF_BUILD_EXAMPLES)

#############################################################################
# Installation
# This code has to be duplicated for install tree
SET(XDMF_INCLUDE_DIRS_CONFIG 
  ${Xdmf_SOURCE_DIR}/libsrc
  ${Xdmf_BINARY_DIR}/libsrc
  )
SET(XDMF_LIBRARY_DIRS_CONFIG 
  ${LIBRARY_OUTPUT_PATH})

IF(XDMF_BUILD_VTK)
  SET(XDMF_INCLUDE_DIRS_CONFIG ${XDMF_INCLUDE_DIRS_CONFIG} 
    ${Xdmf_SOURCE_DIR}/vtk
    ${Xdmf_BINARY_DIR}/vtk
    )
ENDIF(XDMF_BUILD_VTK)

IF(XDMF_SYSTEM_HDF5)
  SET(XDMF_INCLUDE_DIRS_CONFIG ${XDMF_INCLUDE_DIRS_CONFIG}
    ${HDF5_ZLIB_INCLUDE_DIRS})
ELSE(XDMF_SYSTEM_HDF5)
  SET(XDMF_INCLUDE_DIRS_CONFIG ${XDMF_INCLUDE_DIRS_CONFIG}
    ${Xdmf_SOURCE_DIR}/Utilities/hdf5
    ${Xdmf_BINARY_DIR}/Utilities/hdf5)
ENDIF(XDMF_SYSTEM_HDF5)

# Save the compiler settings so another project can import them.
INCLUDE(${CMAKE_ROOT}/Modules/CMakeExportBuildSettings.cmake)
CMAKE_EXPORT_BUILD_SETTINGS(${Xdmf_BINARY_DIR}/XDMFBuildSettings.cmake)

SET(XDMF_BUILD_SETTINGS_FILE ${Xdmf_BINARY_DIR}/XDMFBuildSettings.cmake)
SET(XDMF_USE_FILE "${Xdmf_SOURCE_DIR}/UseXDMF.cmake")
SET(XDMF_KITS_DIR "${Xdmf_BINARY_DIR}/vtk/Utilities")

SET(XDMF_HDF5_LIBRARIES ${HDF5_LIBRARY})

# Save library dependencies.
EXPORT_LIBRARY_DEPENDENCIES(${Xdmf_BINARY_DIR}/XDMFLibraryDepends.cmake)
SET(XDMF_LIBRARY_DEPENDS_FILE ${Xdmf_BINARY_DIR}/XDMFLibraryDepends.cmake)

CONFIGURE_FILE(${Xdmf_SOURCE_DIR}/XDMFConfig.cmake.in
  ${Xdmf_BINARY_DIR}/XDMFConfig.cmake @ONLY IMMEDIATE)

SET(XDMF_INSTALL_BIN_DIR "/lib")
SET(XDMF_INSTALL_LIB_DIR "/lib")
SET(XDMF_INSTALL_INCLUDE_DIR "/include")
SET(XDMF_INSTALL_INCLUDE_VTK_DIR "/include")

# Because INSTALL_* commands require a leading / and because INSTALL (cmake 2.4
# and newer) requires no leading / to install under INSTALL_PREFIX, we
# are stripping the leading /. In the future, there should be no leading
# / in any install directory variables
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_LIB_DIR_CM24 "${XDMF_INSTALL_LIB_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_BIN_DIR_CM24 "${XDMF_INSTALL_BIN_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_INCLUDE_DIR_CM24 "${XDMF_INSTALL_INCLUDE_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_PACKAGE_DIR_CM24 "${XDMF_LIB_INCLUDE_DIR}/vtk")

IF(NOT XDMF_INSTALL_NO_DEVELOPMENT)
    SET(XDMF_INSTALL_LIB_DIR_CM24 ${VTK_INSTALL_LIB_DIR_CM24})
    SET(XDMF_INSTALL_BIN_DIR_CM24 ${VTK_INSTALL_BIN_DIR_CM24})
    SET(XDMF_INSTALL_INCLUDE_DIR_CM24 ${VTK_INSTALL_INCLUDE_DIR_CM24})
    # Re-make XDMFConfig.cmake for installation
   SET(DOLLAR "$")
   SET(XDMF_INSTALL_DIR_CONFIG ${CMAKE_INSTALL_PREFIX})
   SET(XDMF_INCLUDE_DIRS_CONFIG 
     ${DOLLAR}{XDMF_INSTALL_DIR}/include
     )
   SET(XDMF_LIBRARY_DIRS_CONFIG 
     ${DOLLAR}{XDMF_INSTALL_DIR}/lib
     )
   SET(XDMF_USE_FILE ${DOLLAR}{XDMF_INSTALL_DIR}/lib/XdmfCMake/UseXDMF.cmake)
   SET(XDMF_BUILD_SETTINGS_FILE ${DOLLAR}{XDMF_INSTALL_DIR}/lib/XdmfCMake/XDMFBuildSettings.cmake)
   SET(XDMF_LIBRARY_DEPENDS_FILE ${DOLLAR}{XDMF_INSTALL_DIR}/lib/XdmfCMake/XDMFLibraryDepends.cmake)
   SET(XDMF_KITS_DIR ${DOLLAR}{XDMF_INSTALL_DIR}/vtk)
   CONFIGURE_FILE(${Xdmf_SOURCE_DIR}/XDMFConfig.cmake.in
       ${Xdmf_BINARY_DIR}/XDMFConfigInstall.cmake @ONLY IMMEDIATE)
   INSTALL(FILES
       ${Xdmf_BINARY_DIR}/XDMFLibraryDepends.cmake
       ${Xdmf_BINARY_DIR}/XDMFBuildSettings.cmake
       DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24}/XdmfCMake
       COMPONENT Development
       )
   INSTALL(FILES
       ${Xdmf_BINARY_DIR}/XDMFConfigInstall.cmake
       DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24}/XdmfCMake
       COMPONENT Development
       RENAME XDMFConfig.cmake
       )
ENDIF(NOT XDMF_INSTALL_NO_DEVELOPMENT)

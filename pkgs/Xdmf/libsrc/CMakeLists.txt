#-----------------------------------------------------------------------------
# Does Xdmf require support for 64 bit file systems
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CheckCXXSourceRuns.cmake)
FILE(READ "${CMAKE_CURRENT_SOURCE_DIR}/vtkRequireLargeFilesSupport.cxx"
  XDMF_REQUIRE_LARGE_FILE_SUPPORT_FILE)
CHECK_CXX_SOURCE_RUNS("${XDMF_REQUIRE_LARGE_FILE_SUPPORT_FILE}"
  CMAKE_REQUIRE_LARGE_FILE_SUPPORT "Xdmf support for 64 bit file systems")
SET(XDMF_REQUIRE_LARGE_FILE_SUPPORT ${CMAKE_REQUIRE_LARGE_FILE_SUPPORT})
# This Needs to be done since Xdmf/CMakeLists.txt != PV/Utils/Xdmf2/CMakeLists.txt
# Because INSTALL_* commands require a leading / and because INSTALL (cmake 2.4
# and newer) requires no leading / to install under INSTALL_PREFIX, we
# are stripping the leading /. In the future, there should be no leading
# / in any install directory variables
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_LIB_DIR_CM24 "${XDMF_INSTALL_LIB_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_BIN_DIR_CM24 "${XDMF_INSTALL_BIN_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_INCLUDE_DIR_CM24 "${XDMF_INSTALL_INCLUDE_DIR}")
STRING(REGEX REPLACE "^/" "" VTK_INSTALL_PACKAGE_DIR_CM24 "${XDMF_LIB_INCLUDE_DIR}/vtk")
#
SET(XDMF_INSTALL_LIB_DIR_CM24 ${VTK_INSTALL_LIB_DIR_CM24})
SET(XDMF_INSTALL_BIN_DIR_CM24 ${VTK_INSTALL_BIN_DIR_CM24})
SET(XDMF_INSTALL_INCLUDE_DIR_CM24 ${VTK_INSTALL_INCLUDE_DIR_CM24})
SET(XDMF_INSTALL_PACKAGE_DIR_CM24 ${VTK_INSTALL_PACKAGE_DIR_CM24})



SET( XdmfFormatSource
  XdmfArray
  XdmfDOM
  XdmfDataDesc
  XdmfExpression
  XdmfValues
  XdmfValuesXML
  XdmfValuesHDF
  XdmfHeavyData
  XdmfHDF
  XdmfHDFSupport
  XdmfLightData
  XdmfInformation
  XdmfDataItem
  XdmfDataStructure
  XdmfDataTransform
  XdmfElement
  XdmfObject
  )

SET( XdmfExtraSource
  XdmfExprLex
  XdmfExprYacc
  )

SET( XdmfModelSource
  XdmfRoot
  XdmfDomain
  XdmfAttribute
  XdmfGrid
  XdmfTopology
  XdmfTime
  XdmfRegion
  XdmfSet
  XdmfMap
  XdmfGeometry
  XdmfH5Driver
  XdmfDsm
  XdmfDsmBuffer
  XdmfDsmComm
  XdmfDsmMsg
  )

# IF(XDMF_BUILD_MPI)
#     SET(XdmfModelSource ${XdmfModelSource} XdmfH5Driver)
# ENDIF(XDMF_BUILD_MPI) 

SET( XdmfMpiSource
  XdmfDsmCommMpi
  )

SET(XdmfFortranSource
  XdmfFortran
  )

SET( XdmfSTLConverterSource
    XdmfSTLConverter
    )

SET( XdmfMySQLSource
    XdmfValuesMySQL
    )

SET( XdmfPythonSource
  XdmfPython
  )

SET(XDMF_PYTHON_FILE ${Xdmf_SOURCE_DIR}/libsrc/Xdmf.py)
OPTION(XDMF_REGENERATE_WRAPPERS
  "Regenerate XdmfPython.cxx" OFF)
MARK_AS_ADVANCED(XDMF_REGENERATE_WRAPPERS)
IF(XDMF_REGENERATE_WRAPPERS)
  FIND_PROGRAM(SWIG_EXECUTABLE 
    NAMES swig-1.3 swig
    PATHS /usr/bin
    DOC "Swig executable for generating wrapped files")
  MARK_AS_ADVANCED(SWIG_EXECUTABLE)
  IF(SWIG_EXECUTABLE)
    SET(XDMF_SWIG_INCLUDES 
      -I${Xdmf_SOURCE_DIR} 
      -I${Xdmf_SOURCE_DIR}/libsrc 
      -I${Xdmf_BINARY_DIR}/libsrc 
      -I${Xdmf_SOURCE_DIR}/Ice/libsrc
      -I${Xdmf_BINARY_DIR}/Ice/libsrc
      )
    GET_FILENAME_COMPONENT(SWIG_EXE_PATH ${SWIG_EXECUTABLE}
      PATH)
    SET(SWIG_LIB_PATH "${SWIG_EXE_PATH}/Lib")
    # MESSAGE("SWIG_LIB_PATH ${SWIG_LIB_PATH}")
    IF(EXISTS "${SWIG_LIB_PATH}")
      SET(XDMF_SWIG_INCLUDES
        ${XDMF_SWIG_INCLUDES}
        -I${SWIG_LIB_PATH}
        -I${SWIG_LIB_PATH}/tcl
        -I${SWIG_LIB_PATH}/python
        )
    ENDIF(EXISTS "${SWIG_LIB_PATH}")
    FOREACH(LANGUAGE Python)
      STRING(TOLOWER ${LANGUAGE} lang)
      SET(swigtarget ${Xdmf_BINARY_DIR}/libsrc/Xdmf${LANGUAGE}.cxx)
      IF(XDMF_BUILD_MPI)
        ADD_CUSTOM_COMMAND(
         OUTPUT ${swigtarget}
         DEPENDS ${Xdmf_SOURCE_DIR}/libsrc/Xdmf.i
         COMMAND ${SWIG_EXECUTABLE}
         ARGS -v -c++ -make_default -includeall "-${lang}" ${swig_extra_args} ${XDMF_SWIG_INCLUDES} -o "${swigtarget}" Xdmf.i)
      ELSE(XDMF_BUILD_MPI) 
        ADD_CUSTOM_COMMAND(
         OUTPUT ${swigtarget}
         DEPENDS ${Xdmf_SOURCE_DIR}/libsrc/Xdmf.i
         COMMAND ${SWIG_EXECUTABLE}
         ARGS -v -c++ -DXDMF_NO_MPI -make_default -includeall "-${lang}" ${swig_extra_args} ${XDMF_SWIG_INCLUDES} -o "${swigtarget}" Xdmf.i)
      ENDIF(XDMF_BUILD_MPI) 

    ENDFOREACH(LANGUAGE)
    ADD_CUSTOM_TARGET(RerunSwig DEPENDS 
      ${Xdmf_BINARY_DIR}/libsrc/XdmfPython.cxx
      )
  ENDIF(SWIG_EXECUTABLE)
  SET(XDMF_PYTHON_FILE ${Xdmf_BINARY_DIR}/libsrc/Xdmf.py)
ENDIF(XDMF_REGENERATE_WRAPPERS)


OPTION(XDMF_REGENERATE_YACCLEX
  "Regenerate XdmfExprYacc.cxx and XdmfExprLex.cxx" OFF)
MARK_AS_ADVANCED(XDMF_REGENERATE_YACCLEX)
IF(XDMF_REGENERATE_YACCLEX)
  FIND_PROGRAM(YACC_EXECUTABLE
    NAMES yacc bison
    PATHS /usr/bin
    DOC "Yacc or Bison executable")
  FIND_PROGRAM(FLEX_EXECUTABLE
    NAMES flex
    PATHS /usr/bin
    DOC "Flex executable")
  MARK_AS_ADVANCED(YACC_EXECUTABLE FLEX_EXECUTABLE)
  IF(YACC_EXECUTABLE)
    SET(BISON_FLAGS)
    IF(YACC_EXECUTABLE MATCHES "bison")
      SET(BISON_FLAGS "-y")
    ENDIF(YACC_EXECUTABLE MATCHES "bison")
    ADD_CUSTOM_COMMAND(
      OUTPUT ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.tab.c
      DEPENDS ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.y
      COMMAND ${YACC_EXECUTABLE}
      ARGS ${BISON_FLAGS} -d -b ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc
      -p dice_yy
      ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.y
      )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.tab.h
      DEPENDS ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.y
      COMMAND ${YACC_EXECUTABLE}
      ARGS ${BISON_FLAGS} -d -b ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc
      -p dice_yy
      ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.y
      )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.cxx
      DEPENDS ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.tab.c
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.tab.c
      ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.cxx
      )
    ADD_CUSTOM_TARGET(RerunYacc DEPENDS
      ${Xdmf_BINARY_DIR}/libsrc/XdmfExprYacc.cxx)
  ENDIF(YACC_EXECUTABLE)
  IF(FLEX_EXECUTABLE)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${Xdmf_BINARY_DIR}/libsrc/XdmfExprLex.cxx
      DEPENDS ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.l
      COMMAND ${FLEX_EXECUTABLE}
      ARGS -Pdice_yy -o${Xdmf_BINARY_DIR}/libsrc/XdmfExprLex.cxx 
      ${Xdmf_SOURCE_DIR}/libsrc/XdmfExpr.l)
    ADD_CUSTOM_TARGET(RerunLex DEPENDS
      ${Xdmf_BINARY_DIR}/libsrc/XdmfExprLex.cxx)
  ENDIF(FLEX_EXECUTABLE)

ENDIF(XDMF_REGENERATE_YACCLEX)


IF (VTK_USE_MPI)
  # Needed for mpich 2
  ADD_DEFINITIONS("-DMPICH_IGNORE_CXX_SEEK")
ENDIF (VTK_USE_MPI)

INCLUDE (${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
IF(WIN32)
  CHECK_TYPE_SIZE(__int64        SIZEOF___INT64)
ENDIF(WIN32)
CHECK_TYPE_SIZE("long long"    SIZEOF_LONG_LONG)

INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
CHECK_INCLUDE_FILE(malloc.h XDMF_HAVE_MALLOC_H)

# System to Build
STRING(REGEX REPLACE "-" "_" XDMF_SYSTEM ${CMAKE_SYSTEM_NAME})

SET(ARCH_TO_BUILD ${XDMF_SYSTEM} CACHE INTERNAL "Host Arcitecture : Linux IRIXN32 IRIX64 AIX CYGWIN")
# Allow the user to customize their build with some local options
#
ADD_DEFINITIONS(-D${ARCH_TO_BUILD})

SET(XDMF_SIZEOF_CHAR   ${CMAKE_SIZEOF_CHAR})
SET(XDMF_SIZEOF_DOUBLE ${CMAKE_SIZEOF_DOUBLE})
SET(XDMF_SIZEOF_FLOAT  ${CMAKE_SIZEOF_FLOAT})
SET(XDMF_SIZEOF_INT    ${CMAKE_SIZEOF_INT})
SET(XDMF_SIZEOF_LONG   ${CMAKE_SIZEOF_LONG})
SET(XDMF_SIZEOF_SHORT  ${CMAKE_SIZEOF_SHORT})
INCLUDE (${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
IF(WIN32)
  CHECK_TYPE_SIZE(__int64        SIZEOF___INT64)
ENDIF(WIN32)
CHECK_TYPE_SIZE("long long"    SIZEOF_LONG_LONG)
CHECK_FUNCTION_EXISTS(strtoll HAVE_STRTOLL)

# Don't Really Need to check these ...
SET(XDMF_VOID_VALID 1)
SET(XDMF_VOID "void")
SET(XDMF_PTR_VALID 1)
SET(XDMF_PTR "void *")
SET(XDMF_CHAR_VALID 1)
SET(XDMF_CHAR "char")
SET(XDMF_8_INT_VALID 1)
SET(XDMF_8_INT  "char")
SET(XDMF_8_U_INT  "unsigned char")
SET(XDMF_16_INT_VALID 1)
SET(XDMF_16_INT "short")
SET(XDMF_16_U_INT "unsigned short")
# These should be Valid
IF(${XDMF_SIZEOF_FLOAT} MATCHES 4)
  SET(XDMF_32_FLOAT_VALID 1)
  SET(XDMF_FLOAT "float")
  SET(XDMF_32_FLOAT "float")
ELSE (${XDMF_SIZEOF_FLOAT} MATCHES 4)
  MESSAGE(SEND_ERROR "Can't Find a 32 Bit Float")
ENDIF (${XDMF_SIZEOF_FLOAT} MATCHES 4)
IF(${XDMF_SIZEOF_DOUBLE} MATCHES 8)
  SET(XDMF_64_FLOAT_VALID 1)
  SET(XDMF_DOUBLE "double")
  SET(XDMF_64_FLOAT "double")
ELSE (${XDMF_SIZEOF_DOUBLE} MATCHES 8)
  MESSAGE(SEND_ERROR "Can't Find a 64 Bit Float")
ENDIF (${XDMF_SIZEOF_DOUBLE} MATCHES 8)
# These are sometimes different
IF(${XDMF_SIZEOF_INT} MATCHES 4)
  SET(XDMF_32_INT_VALID 1)
  SET(XDMF_32_INT "int")
  SET(XDMF_32_U_INT "unsigned int")
  SET(XDMF_32_S_INT "int")
ELSE (${XDMF_SIZEOF_INT} MATCHES 4)
  IF(${XDMF_SIZEOF_LONG} MATCHES 4)
    SET(XDMF_32_INT_VALID 1)
    SET(XDMF_32_INT "long")
    SET(XDMF_32_U_INT "unsigned long")
    SET(XDMF_32_S_INT "long")
  ELSE(${XDMF_SIZEOF_LONG} MATCHES 4)
    MESSAGE(SEND_ERROR "Can't Find a 32 Bit Integer")
  ENDIF(${XDMF_SIZEOF_LONG} MATCHES 4)
ENDIF(${XDMF_SIZEOF_INT} MATCHES 4)

INCLUDE(${Xdmf_SOURCE_DIR}/CMake/CheckFor64BitStreams.cmake)
CHECK_FOR_64BIT_STREAMS(XDMF_HAVE_64BIT_STREAMS)

# Find include files
FIND_PATH(XDMF_HAVE_FCNTL fcntl.h /usr/include /usr/include/sys)
FIND_PATH(XDMF_HAVE_NETINET in.h /usr/include/netinet /usr/include /usr/include/sys)
FIND_PATH(XDMF_HAVE_MMAN mman.h /usr/include/sys /usr/include)
MARK_AS_ADVANCED(XDMF_HAVE_FCNTL XDMF_HAVE_NETINET XDMF_HAVE_MMAN)


##
CONFIGURE_FILE(${Xdmf_SOURCE_DIR}/libsrc/XdmfConfig.h.in
  ${Xdmf_BINARY_DIR}/libsrc/XdmfConfig.h)

IF(NOT WIN32)
  ADD_DEFINITIONS(-D_HPUX_SOURCE)
ENDIF(NOT WIN32)

SET(XdmfSources ${XdmfFormatSource} ${XdmfExtraSource} ${XdmfModelSource})

OPTION(XDMF_USE_MYSQL
    "Build Support for MySQL DataItems" OFF)
MARK_AS_ADVANCED(XDMF_USE_MYSQL)
IF(XDMF_USE_MYSQL)
    FIND_PATH(MYSQL_INCLUDE_PATH
        mysql.h
        /usr/include
        /usr/include/mysql
        /usr/local/include
        /usr/local/include/mysql)
    FIND_LIBRARY(MYSQL_CLIENT_LIBRARY
        mysqlclient
        /usr/lib
        /usr/lib/mysql
        /usr/local/lib
        /usr/local/lib/mysql)
    SET(XdmfSources ${XdmfSources} ${XdmfMySQLSource})
    ADD_DEFINITIONS("-DXDMF_USE_MYSQL")
    INCLUDE_DIRECTORIES(${MYSQL_INCLUDE_PATH})
ENDIF(XDMF_USE_MYSQL)

OPTION(XDMF_USE_STLCONVERTER "Build Special STL Converter" OFF)
MARK_AS_ADVANCED(XDMF_USE_STLCONVERTER)
IF(XDMF_USE_STLCONVERTER)
  INSTALL(
      FILES ${Xdmf_SOURCE_DIR}/libsrc/XdmfSTLConverter.txx
      DESTINATION ${XDMF_INSTALL_INCLUDE_DIR_CM24})
ENDIF(XDMF_USE_STLCONVERTER)


## MPI ###
IF(XDMF_BUILD_MPI)
    SET(XdmfSources ${XdmfSources} ${XdmfMpiSource})
ELSE(XDMF_BUILD_MPI)
    ADD_DEFINITIONS("-DXDMF_NO_MPI")
ENDIF(XDMF_BUILD_MPI)

IF(VTK_BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)
  SET(LIBTYPE SHARED)
ELSE(VTK_BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)
  SET(LIBTYPE STATIC)
ENDIF(VTK_BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)

ADD_LIBRARY(Xdmf ${LIBTYPE} ${XdmfSources})
ADD_LIBRARY(XdmfFortran ${LIBTYPE} ${XdmfFortranSource})

IF(XDMF_USE_MYSQL)
    ADD_DEFINITIONS("-DXDMF_USE_MYSQL")
    INCLUDE_DIRECTORIES(${MYSQL_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(Xdmf ${MYSQL_CLIENT_LIBRARY} )
ENDIF(XDMF_USE_MYSQL)

IF(XDMF_BUILD_MPI)
  IF(MPI_EXTRA_LIBRARY)
    TARGET_LINK_LIBRARIES(Xdmf ${MPI_EXTRA_LIBRARY})
  ENDIF(MPI_EXTRA_LIBRARY)
  TARGET_LINK_LIBRARIES(Xdmf ${MPI_LIBRARY} )
ENDIF(XDMF_BUILD_MPI)

TARGET_LINK_LIBRARIES(Xdmf 
  ${HDF5_LIBRARY} 
  ${XDMF_MPI_LIBRARIES} 
  ${LIBXML2_LIBRARY} 
  ${XDMF_ZLIB_LIBRARIES})

TARGET_LINK_LIBRARIES(XdmfFortran Xdmf)

OPTION(XDMF_WRAP_PYTHON "Build a Python wrapping interface" OFF)
SET (temp_VTK_PYTHON FALSE)
IF (VTK_WRAP_PYTHON AND NOT PARAVIEW_ENABLE_PYTHON)
  SET (temp_VTK_PYTHON TRUE)
ENDIF (VTK_WRAP_PYTHON AND NOT PARAVIEW_ENABLE_PYTHON)

IF(XDMF_WRAP_PYTHON OR temp_VTK_PYTHON)
  FIND_PACKAGE(PythonExecutable)
  FIND_PACKAGE(PythonLibs)
ENDIF(XDMF_WRAP_PYTHON OR temp_VTK_PYTHON)

IF(XDMF_WRAP_PYTHON)
  IF(NOT BUILD_SHARED_LIBS)
    MESSAGE("Python wrapping only works on shared libraries")
  ENDIF(NOT BUILD_SHARED_LIBS)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  ADD_LIBRARY(_Xdmf SHARED ${XdmfPythonSource})
  TARGET_LINK_LIBRARIES(_Xdmf Xdmf ${HDF5_LIBRARY} ${XDMF_ZLIB_LIBRARIES})
  IF(WIN32)
    TARGET_LINK_LIBRARIES(_Xdmf ${PYTHON_LIBRARY})
  ENDIF(WIN32)
  SET(PYTHON_VERSION_NAME "python")
  GET_FILENAME_COMPONENT(PYTHON_VERSION_NAME ${PYTHON_INCLUDE_PATH} NAME)
  SET_TARGET_PROPERTIES(_Xdmf PROPERTIES PREFIX "")
  ADD_DEPENDENCIES(_Xdmf Xdmf)
  INSTALL(FILES ${Xdmf_SOURCE_DIR}/libsrc/__init__.py 
      DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24}/${PYTHON_VERSION_NAME}/site-packages/Xdmf
      COMPONENT Runtime
      )
    INSTALL(FILES ${XDMF_PYTHON_FILE}
        DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24}/${PYTHON_VERSION_NAME}/site-packages/Xdmf
        COMPONENT Runtime
        )
ENDIF(XDMF_WRAP_PYTHON)

IF(XDMF_BUILD_TESTING)
  SUBDIRS(Testing)
ENDIF(XDMF_BUILD_TESTING)

IF(NOT XDMF_INSTALL_NO_DEVELOPMENT)
  FILE(GLOB devFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
  INSTALL(
      FILES ${devFiles}
      DESTINATION ${XDMF_INSTALL_INCLUDE_DIR_CM24})
  INSTALL(
      FILES ${Xdmf_BINARY_DIR}/libsrc/XdmfConfig.h
      DESTINATION ${XDMF_INSTALL_INCLUDE_DIR_CM24})
ENDIF(NOT XDMF_INSTALL_NO_DEVELOPMENT)

IF(NOT XDMF_INSTALL_NO_LIBRARIES)
  INSTALL(TARGETS Xdmf XdmfFortran
    RUNTIME DESTINATION ${XDMF_INSTALL_BIN_DIR_CM24} COMPONENT Runtime
    LIBRARY DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24} COMPONENT Runtime
    ARCHIVE DESTINATION ${XDMF_INSTALL_LIB_DIR_CM24} COMPONENT Development)
  IF(XDMF_WRAP_PYTHON)
    INSTALL(
        TARGETS _Xdmf
        DESTINATION
        ${XDMF_INSTALL_LIB_DIR_CM24}/${PYTHON_VERSION_NAME}/lib-dynload)
  ENDIF(XDMF_WRAP_PYTHON)
ENDIF(NOT XDMF_INSTALL_NO_LIBRARIES)

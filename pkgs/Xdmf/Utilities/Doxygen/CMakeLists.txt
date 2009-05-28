#
# Build the documentation
#

INCLUDE (${CMAKE_ROOT}/Modules/Documentation.cmake)

IF (BUILD_DOCUMENTATION)
  OPTION (DOCUMENTATION_DOWNLOAD_VTK_TAGFILE
         "Download the VTK tag file." ON)
  MARK_AS_ADVANCED (
    DOCUMENTATION_DOWNLOAD_VTK_TAGFILE
  )


  #
  # Configure the script and the doxyfile, then add target
  #
  CONFIGURE_FILE(
    ${Xdmf_SOURCE_DIR}/Utilities/Doxygen/doxyfile.in
    ${Xdmf_BINARY_DIR}/Utilities/Doxygen/doxyfile)

  CONFIGURE_FILE(
    ${Xdmf_SOURCE_DIR}/Utilities/Doxygen/doc_makeall.sh.in
    ${Xdmf_BINARY_DIR}/Utilities/Doxygen/doc_makeall.sh)

  ADD_CUSTOM_TARGET(XdmfDoxygenDoc
    ${BASH}
    ${Xdmf_BINARY_DIR}/Utilities/Doxygen/doc_makeall.sh)

ENDIF (BUILD_DOCUMENTATION)

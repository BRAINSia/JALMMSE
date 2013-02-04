include(ExternalProject)

#
# argv1 == program name
# argv2 == src variable name
# argv2 == xml name
macro(STANDALONEGENERATEMODULESCRIPT progname)
  #
  # Next, do the configure thing.
  #
  # get the directory to throw the dummy into
  # message(STATUS "PROGRAM NAME ${progname}")
  get_target_property(BGCbin_dir ${progname} RUNTIME_OUTPUT_DIRECTORY)

  set(BGCbin_name ${progname})
  IF(APPLE)
    SET(OS_VARNAME_FOR_LIBRARY_PATH "DYLD_LIBRARY_PATH")
  ELSE(APPLE)
    SET(OS_VARNAME_FOR_LIBRARY_PATH "LD_LIBRARY_PATH")
  ENDIF(APPLE)

  CONFIGURE_FILE(${STANDALONE_CMAKE_HELPER_DIR}/Module_Dummy.in ${BGCbin_dir}/Modules/${progname} )

# FILE(WRITE ${BGCbin_dir}/Modules/${progname}
# "#!/usr/bin/tclsh
# catch {set script [info script]}
# catch {set script [file normalize \$script]}
# catch {set execdir [file dirname [file dirname \$script ]]}
# set env(${OS_VARNAME_FOR_LIBRARY_PATH}) [ exec \$execdir/brains3_setup.sh ${OS_VARNAME_FOR_LIBRARY_PATH} ]
# set command \"\$execdir/${BGCbin_name} \$argv\"
# set fp [ open \"| \$command |& cat\" \"r\"]
# while { ![eof \$fp ] } {
#     gets \$fp line
#     puts \$line
# }
# if { [catch \"close \$fp\" res] } {
#     exit [ lindex \$errorCode 2 ]
# } else {
#     exit 0
# }
# "
# )


  install(PROGRAMS ${BGCbin_dir}/Modules/${progname}
    DESTINATION bin/Modules PERMISSIONS WORLD_EXECUTE)
endmacro(STANDALONEGENERATEMODULESCRIPT)


#-----------------------------------------------------------------------------
# Build the optional DEBUGIMAGEVIEWER
if(NOT SETOPTIONALDEBUGIMAGEVIEWER)
macro(SETOPTIONALDEBUGIMAGEVIEWER)
if(STANDALONE_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" ON)
else(STANDALONE_BUILD)
  option(USE_DEBUG_IMAGE_VIEWER "Use the DEBUG_IMAGE_VIEWER for debugging" OFF)
endif(STANDALONE_BUILD)

mark_as_advanced(USE_DEBUG_IMAGE_VIEWER)
set(OPTIONAL_DEBUG_LINK_LIBRARIES) ## Set it to empty as the default
if( USE_DEBUG_IMAGE_VIEWER )
   if(NOT KWWidgets_SOURCE_DIR)
     find_package(KWWidgets REQUIRED)
     include(${KWWidgets_USE_FILE})
   endif(NOT KWWidgets_SOURCE_DIR)
   add_definitions(-DUSE_DEBUG_IMAGE_VIEWER)
   find_path(DEBUG_IMAGE_VIEWER_INCLUDE_DIR DebugImageViewerClient.h ${CMAKE_INSTALL_PREFIX}/include)
   include_directories(${DEBUG_IMAGE_VIEWER_INCLUDE_DIR})
   set(OPTIONAL_DEBUG_LINK_LIBRARIES ${KWWidgets_LIBRARIES})
endif( USE_DEBUG_IMAGE_VIEWER )
endmacro(SETOPTIONALDEBUGIMAGEVIEWER)
endif(NOT SETOPTIONALDEBUGIMAGEVIEWER)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create executables for Slicer or STANDALONE3
if(NOT CONFIGURESTANDALONEORSLICERPROPERTIES)
macro(CONFIGURESTANDALONEORSLICERPROPERTIES PROGNAME PROGCLI PROGSOURCES EXTRA_LIBS EXTRA_HEADERS)

  find_package(GenerateCLP NO_MODULE REQUIRED)
  include(${GenerateCLP_USE_FILE})

  get_filename_component(TMP_FILENAME ${PROGCLI} NAME_WE)
  set(PROGCLI_HEADER "${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h")

  set(CLP_SOURCES ${PROGSOURCES})
  if(EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/ModuleLogo.h)
    GENERATECLP(CLP_SOURCES ${PROGCLI} ${CMAKE_CURRENT_SOURCE_DIR}/ModuleLogo.h)
  else()
    GENERATECLP(CLP_SOURCES ${PROGCLI} )
  endif()
  add_executable( ${PROGNAME} ${CLP_SOURCES} ${EXTRA_HEADERS} ${PROGCLI_HEADER})
  if(WIN32)
    set(STANDALONE_ITK_LIBS "")
  else(WIN32)
    set(STANDALONE_ITK_LIBS ${ITK_LIBRARIES})
  endif(WIN32)
  target_link_libraries (${PROGNAME} ${STANDALONE_ITK_LIBS} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

  if (Slicer_SOURCE_DIR)
    ### If building as part of the Slicer_SOURCE_DIR, then only build the shared object, and not the command line program.

    add_library(${PROGNAME}Lib SHARED ${CLP_SOURCES} ${PROGCLI_HEADER})
    set_target_properties (${PROGNAME}Lib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
    slicer3_set_plugins_output_path(${PROGNAME}Lib)
    target_link_libraries (${PROGNAME}Lib ${STANDALONE_ITK_LIBS} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

    # install each target in the production area (where it would appear in an
    # installation) and install each target in the developer area (for running
    # from a build)
    slicer3_set_plugins_output_path(${PROGNAME})
    set(TARGETS ${PROGNAME}Lib ${PROGNAME})
    slicer3_install_plugins(${TARGETS})
  else (Slicer_SOURCE_DIR)
    ### If building outside of Slicer3, then only build the command line executable.
    IF(STANDALONE_BUILD)
      STANDALONEGENERATEMODULESCRIPT(${PROGNAME})
    ENDIF(STANDALONE_BUILD)
    INSTALL(TARGETS ${PROGNAME}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
  endif (Slicer_SOURCE_DIR)

endmacro(CONFIGURESTANDALONEORSLICERPROPERTIES PROGNAME PROGCLI PROGSOURCES)
endif(NOT CONFIGURESTANDALONEORSLICERPROPERTIES)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
## A macro to create CLP dependant libraries for Slicer or STANDALONE3
if(NOT CONFIGURESTANDALONEORSLICERLIBRARY)
macro(CONFIGURESTANDALONEORSLICERLIBRARY LIBNAME LIBCLI LIBSOURCES EXTRA_LIBS)

  find_package(GenerateCLP NO_MODULE REQUIRED)
  include(${GenerateCLP_USE_FILE})

  get_filename_component(TMP_FILENAME ${LIBCLI} NAME_WE)
  set(LIBCLI_HEADER "${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h")

  set(CLP_SOURCES ${LIBSOURCES})
  if(EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/ModuleLogo.h)
    GENERATECLP(CLP_SOURCES ${LIBCLI} ${CMAKE_CURRENT_SOURCE_DIR}/ModuleLogo.h)
  else()
    GENERATECLP(CLP_SOURCES ${LIBCLI} )
  endif()
  add_library( ${LIBNAME} ${CLP_SOURCES} ${LIBCLI_HEADER})
  target_link_libraries (${LIBNAME} ${ITK_LIBRARIES} ${OPTIONAL_DEBUG_LINK_LIBRARIES} ${EXTRA_LIBS} )

  if (Slicer_SOURCE_DIR)
    ### If building as part of the Slicer_SOURCE_DIR, then only build the shared object, and not the command line program.
    # install each target in the production area (where it would appear in an
    # installation) and install each target in the developer area (for running
    # from a build)
    slicer3_set_plugins_output_path(${LIBNAME})
    set(TARGETS ${LIBNAME})
    slicer3_install_plugins(${TARGETS})
  else (Slicer_SOURCE_DIR)
    ### If building outside of Slicer3, then only build the command line executable.
    IF(STANDALONE_BUILD)
      STANDALONEGENERATEMODULESCRIPT(${LIBNAME})
    ENDIF(STANDALONE_BUILD)
    INSTALL(TARGETS ${LIBNAME}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
  endif (Slicer_SOURCE_DIR)

endmacro(CONFIGURESTANDALONEORSLICERLIBRARY LIBNAME LIBCLI LIBSOURCES)
endif(NOT CONFIGURESTANDALONEORSLICERLIBRARY)

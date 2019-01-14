include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------


################################################################################
#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

#-----------------------------------------------------------------------------
find_package(ITK REQUIRED)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

# For Slicer 4 builds, simply call the standard macro:
if(Slicer_BINARY_DIR) # This variable exits only for Slicer4

    SEMMacroBuildCLI(
       NAME ${PROJECT_NAME}
       LOGO_HEADER ${PROJECT_SOURCE_DIR}/ModuleLogo.h
       TARGET_LIBRARIES ${ITK_LIBRARIES} ModuleDescriptionParser
       LINK_DIRECTORIES ${ModuleDescriptionParser_BINARY_DIR}
       INCLUDE_DIRECTORIES ${SlicerBaseCLI_SOURCE_DIR} ${SlicerBaseCLI_BINARY_DIR}
       EXECUTABLE_ONLY
       RUNTIME_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY}
       LIBRARY_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY}
       ARCHIVE_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY}
       INSTALL_RUNTIME_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}
       INSTALL_LIBRARY_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION}
       INSTALL_ARCHIVE_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION}
    )

else(Slicer_BINARY_DIR)
SEMMacroBuildCLI(
       NAME ${PROJECT_NAME}
       EXECUTABLE_ONLY
       LOGO_HEADER ${PROJECT_SOURCE_DIR}/ModuleLogo.h
       TARGET_LIBRARIES ${ITK_LIBRARIES} ModuleDescriptionParser
       LINK_DIRECTORIES ${ModuleDescriptionParser_BINARY_DIR}
       INCLUDE_DIRECTORIES ${SlicerBaseCLI_SOURCE_DIR} ${Slicer_SOURCE_DIR}/Applications/CLI/DiffusionApplications/DiffusionApplicationsCommon
       RUNTIME_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY}
       LIBRARY_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY}
       ARCHIVE_OUTPUT_DIRECTORY ${${CMAKE_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY}
       INSTALL_RUNTIME_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}
       INSTALL_LIBRARY_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION}
       INSTALL_ARCHIVE_DESTINATION ${${CMAKE_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION}
)
endif(Slicer_BINARY_DIR)

if (BUILD_TESTING)
add_subdirectory(TestSuite)
endif (BUILD_TESTING)




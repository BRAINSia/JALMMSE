
add_executable(${CLP}_Compare ${CLP}_Compare.cxx)
target_link_libraries(${CLP}_Compare ${ITK_LIBRARIES} )

set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/TestData)
set(OUTPUT_DIR    ${CMAKE_CURRENT_BINARY_DIR} )
set(LAUNCH_EXE "")

if (Slicer_SOURCE_DIR AND Slicer_BINARY_DIR)
  # And this is for Slicer4. It requires calling the CLI with the full path
  set(Slicer_EXE ${Slicer_BINARY_DIR}/Slicer)
  set(LAUNCH_EXE "${Slicer_EXE} --launch ")
endif ()

# This runs the module with testing parameters:
add_test(NAME ${PROJECT_NAME}_SampleTest1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}>
  --overrideNoise
  --noiseLevel 5
  --ng 0
  --rf 1,1,0
  --re 1,1,0
  --h 999999999
  ${TEST_DATA_DIR}/sampleDWI.nrrd
  ${OUTPUT_DIR}/filtered1DWI.nrrd
)

add_test(NAME ${PROJECT_NAME}_SampleTest2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}>
  --overrideNoise
  --noiseLevel 10
  --ng 4
  --rf 2,2,1
  --re 2,2,1
  --h 2.0
  --maskImage ${TEST_DATA_DIR}/mask.nrrd
  --setZeroBck
  ${TEST_DATA_DIR}/sampleDWI.nrrd
  ${OUTPUT_DIR}/filtered2DWI.nrrd
)

add_test(NAME ${PROJECT_NAME}_SampleTest3
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}>
  --overrideNoise
  --noiseLevel 15
  --ng 15
  --rf 3,3,3
  --re 2,2,1
  --h 1.2
  --onlyUNLM
  ${TEST_DATA_DIR}/sampleDWI.nrrd
  ${OUTPUT_DIR}/filtered3DWI.nrrd
)

add_test(NAME ${PROJECT_NAME}_SampleTest4
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}>
  --ng 1
  --rf 3,3,3
  --re 1,1,0
  --h 2.0
  --onlyUNLM
  --maskImage ${TEST_DATA_DIR}/mask.nrrd
  ${TEST_DATA_DIR}/sampleDWI.nrrd
  ${OUTPUT_DIR}/filtered4DWI.nrrd
)

# This runs a program that checks the output and compares to the one expected:
add_test(NAME ${PROJECT_NAME}_SampleCompare1
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}_Compare>
  ${OUTPUT_DIR}/filtered1DWI.nrrd
  ${TEST_DATA_DIR}/filtered1DWI.nrrd
  1e-9
)
set_property(TEST ${PROJECT_NAME}_SampleCompare1 APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleTest1)

add_test(NAME ${PROJECT_NAME}_SampleCompare2
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}_Compare>
  ${OUTPUT_DIR}/filtered2DWI.nrrd
  ${TEST_DATA_DIR}/filtered2DWI.nrrd
  1e-9
)
set_property(TEST ${PROJECT_NAME}_SampleCompare2 APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleTest2)

add_test(NAME ${PROJECT_NAME}_SampleCompare3
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}_Compare>
  ${OUTPUT_DIR}/filtered3DWI.nrrd
  ${TEST_DATA_DIR}/filtered3DWI.nrrd
  1e-9
)
set_property(TEST ${PROJECT_NAME}_SampleCompare3 APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleTest3)

add_test(NAME ${PROJECT_NAME}_SampleCompare4
  COMMAND ${LAUNCH_EXE} $<TARGET_FILE:${CLP}_Compare>
  ${OUTPUT_DIR}/filtered4DWI.nrrd
  ${TEST_DATA_DIR}/filtered4DWI.nrrd
  1e-9
)
set_property(TEST ${PROJECT_NAME}_SampleCompare4 APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleTest4)

# Clean the outputs created
add_test(NAME ${PROJECT_NAME}_CleanResults
  COMMAND ${CMAKE_COMMAND} -E remove
  ${OUTPUT_DIR}/filtered1DWI.nrrd
  ${OUTPUT_DIR}/filtered2DWI.nrrd
  ${OUTPUT_DIR}/filtered3DWI.nrrd
  ${OUTPUT_DIR}/filtered4DWI.nrrd
)
set_property(TEST ${PROJECT_NAME}_CleanResults APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleCompare1)
set_property(TEST ${PROJECT_NAME}_CleanResults APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleCompare2)
set_property(TEST ${PROJECT_NAME}_CleanResults APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleCompare3)
set_property(TEST ${PROJECT_NAME}_CleanResults APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SampleCompare4)


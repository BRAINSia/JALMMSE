

  add_executable(${CLP}_Compare ${CLP}_Compare.cxx)
  target_link_libraries(${CLP}_Compare ${ITK_LIBRARIES} )

  set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/TestData)
  if (Slicer_SOURCE_DIR)

        # This is a Slicer build
        # Depending on the actual version (3 or 4), the calling sequence differs:
        if(Slicer3_BINARY_DIR)
            # This is for Slicer3
            set(OUTPUT_DIR ${Slicer3_BINARY_DIR}/Testing)
            set(Slicer_EXE ${Slicer3_BINARY_DIR}/Slicer3)
            SET(RUNEXEC
                ${Slicer_EXE}
                --launch
                ${CLP}
            )
            SET(COMPEXEC
                ${Slicer_EXE}
                --launch
                ${CLP}_Compare
            )
        endif(Slicer3_BINARY_DIR)

        if(Slicer_BINARY_DIR)
            # And this is for Slicer4. It requires calling the CLI with the full path
            set(OUTPUT_DIR ${Slicer_BINARY_DIR}/Testing)
            set(Slicer_EXE ${Slicer_BINARY_DIR}/Slicer)
            SET(RUNEXEC
                ${Slicer_EXE}
                --launch
                "${Slicer_BINARY_DIR}/lib/Slicer-${Slicer_VERSION_MAJOR}.${Slicer_VERSION_MINOR}/cli-modules/${CLP}"
            )
            SET(COMPEXEC
                ${Slicer_EXE}
                --launch
                "${Slicer_BINARY_DIR}/bin/${CLP}_Compare"
            )
        endif(Slicer_BINARY_DIR)

  else( Slicer_SOURCE_DIR )
    set( OUTPUT_DIR ${CMAKE_BINARY_DIR}/Testing )
    SET( RUNEXEC ${CMAKE_BINARY_DIR}/${CLP} )
    SET( COMPEXEC ${CMAKE_BINARY_DIR}/${CLP}_Compare )
  endif (Slicer_SOURCE_DIR)

  # This runs the module with testing parameters:
  add_test(${PROJECT_NAME}_SampleTest1
    ${RUNEXEC}
    --overrideNoise
    --noiseLevel 5
    --ng 0
    --rf 1,1,0
    --re 1,1,0
    --h 999999999
    ${TEST_DATA_DIR}/sampleDWI.nrrd
    ${OUTPUT_DIR}/filtered1DWI.nrrd
  )

  add_test(${PROJECT_NAME}_SampleTest2
    ${RUNEXEC}
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

  add_test(${PROJECT_NAME}_SampleTest3
    ${RUNEXEC}
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

  add_test(${PROJECT_NAME}_SampleTest4
    ${RUNEXEC}
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
  add_test(${PROJECT_NAME}_SampleCompare1
    ${COMPEXEC}
    ${OUTPUT_DIR}/filtered1DWI.nrrd
    ${TEST_DATA_DIR}/filtered1DWI.nrrd
    1e-9
  )

  add_test(${PROJECT_NAME}_SampleCompare2
    ${COMPEXEC}
    ${OUTPUT_DIR}/filtered2DWI.nrrd
    ${TEST_DATA_DIR}/filtered2DWI.nrrd
    1e-9
  )

  add_test(${PROJECT_NAME}_SampleCompare3
    ${COMPEXEC}
    ${OUTPUT_DIR}/filtered3DWI.nrrd
    ${TEST_DATA_DIR}/filtered3DWI.nrrd
    1e-9
  )

  add_test(${PROJECT_NAME}_SampleCompare4
    ${COMPEXEC}
    ${OUTPUT_DIR}/filtered4DWI.nrrd
    ${TEST_DATA_DIR}/filtered4DWI.nrrd
    1e-9
  )

  # Clean the outputs created
  add_test(${PROJECT_NAME}_CleanResults
    rm
    ${OUTPUT_DIR}/filtered1DWI.nrrd
    ${OUTPUT_DIR}/filtered2DWI.nrrd
    ${OUTPUT_DIR}/filtered3DWI.nrrd
    ${OUTPUT_DIR}/filtered4DWI.nrrd
  )


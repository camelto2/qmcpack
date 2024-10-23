if (QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  function(ADD_DIRAC_TEST)
  
  endfunction()

  function(ADD_DIRAC_CONVERT_TEST)

  endfunction()

  function(RUN_DIRAC_TEST)

  endfunction()
else(QMC_NO_SLOW_CUSTOM_TESTING_COMMNANDS)
  function(
    ADD_DIRAC_TEST
    TESTNAME
    PROCS
    LAUNCHER
    WORKDIR)
    if (HAVE_MPI)
      add_test(NAME ${TESTNAME} COMMAND ${LAUNCHER} --mpi=${MPIEXEC_NUMPROC_FLAG} --inp=input.inp --mol=system.mol)
    else()
      add_test(NAME ${TESTNAME} COMMAND ${LAUNCHER} --inp=input.inp --mol=system.mol)
    endif()
    set_tests_properties(
      ${TESTNAME}
      PROPERTIES ENVIRONMENT
                 OMP_NUM_THREADS=1
                 PROCESSORS
                 ${PROCS}
                 PROCESSOR_AFFINITY
                 TRUE
                 WORKING_DIRECTORY
                 ${WORKDIR})
    set_property(
      TEST ${TESTNAME}
      APPEND
      PROPERTY LABELS "converter")
  endfunction()

  function(ADD_DIRAC_CONVERT_TEST TESTNAME PREFIX WORKDIR)
    add_test(NAME ${TESTNAME} COMMAND $<TARGET_FILE:convert4qmc> -nojastrow -dirac input_system.out -prefix ${PREFIX})
    set_tests_properties(${TESTNAME} PROPERTIES WORKING_DIRECTORY ${WORKDIR})
    set_property(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;dirac")
  endfunction()

  function(
    RUN_DIRAC_TEST
    BASE_NAME
    SRC_DIR
    NPROCS
    TEST_NAME
  )
    set(FULL_NAME ${BASE_NAME}-np-${NPROCS})
    set(${TEST_NAME} ${FULL_NAME} PARENT_SCOPE)
    set(MY_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME})
    message(VERBOSE "Adding test ${FULL_NAME}")
    copy_directory("${SRC_DIR}" "${MY_WORKDIR}")
    message("workdir: ${MY_WORKDIR}")
    add_dirac_test(${FULL_NAME}-scf NPROCS ${DIRAC_LAUNCHER} ${MY_WORKDIR})
    add_dirac_convert_test(${FULL_NAME}-dirac2qmc ${BASE_NAME} ${MY_WORKDIR})

    set_tests_properties(${FULL_NAME}-dirac2qmc PROPERTIES DEPENDS ${FULL_NAME}-nscf) 

  endfunction()
endif(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)

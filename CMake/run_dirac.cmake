if (QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)
  function(ADD_DIRAC_TEST)
  
  endfunction()
else(QMC_NO_SLOW_CUSTOM_TESTING_COMMNANDS)
  function(
    ADD_DIRAC_TEST
    TESTNAME
    PROCS
    LAUNCHER
    WORKDIR
    TEST_INPUT
    TEST_MOL)
    if (HAVE_MPI)
      add_test(NAME ${TESTNAME} COMMAND ${LAUNCHER} --mpi=${MPIEXEC_NUMPROC_FLAG} --inp=${TEST_INPUT} --mol=${TEST_MOL})
    else()
      add_test(NAME ${TESTNAME} COMMAND ${LAUNCHER} --inp=${TEST_INPUT} --mol=${TEST_MOL})
    endif()
    set_test_properties(
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

endif(QMC_NO_SLOW_CUSTOM_TESTING_COMMANDS)


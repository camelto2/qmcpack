message("Adding DIRAC tests")

if (NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 16)
  message("TEST_MAX_PROCS was unset. Set to 16.")
endif()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()


if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  include("${qmcpack_SOURCE_DIR}/CMake/run_dirac.cmake")
  set(LAST_TEST_NAME "NONE")
  set(DIRAC_TEST_NAME "NONE")
  
  set(DIRAC_SYSTEM "Bi_cosci_ccECP")
  run_dirac_test(
    dirac-${DIRAC_SYSTEM} 
    ${qmcpack_SOURCE_DIR}/examples/molecules/${DIRAC_SYSTEM}
    12 
    DIRAC_TEST_NAME)

endif()

# Locate DIRAC
# Take DIRAC_BIN as hint for location

find_program(DIRAC_LAUNCHER pam HINTS ${DIRAC_BIN})

set(DIRAC_FOUND FALSE)
if(DIRAC_LAUNCHER)
  message(STATUS "DIRAC_LAUNCHER=${DIRAC_EXE}")
  set(DIRAC_FOUND TRUE)
endif()

mark_as_advanced(DIRAC_LAUNCHER DIRAC_FOUND)

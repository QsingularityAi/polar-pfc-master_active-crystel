set(INCLUDE_DIRS "")
set(LIBRARIES "")

# find SHTNS libraries
find_library(_SHTNS_LIB shtns_omp HINTS ${SHTNS_LIB} $ENV{SHTNS_LIB})
if(_SHTNS_LIB)
  get_filename_component(_SHTNS_DIR ${_SHTNS_LIB} PATH CACHE)
  set(LIBRARY_DIRS ${_SHTNS_DIR})
  set(LIBRARIES "${_SHTNS_LIB}")

  # find shtns finclude directory
  find_path(_SHTNS_INC shtns.h HINTS ${SHTNS_INC} $ENV{SHTNS_INC} ${_SHTNS_DIR}/../include)
  if(_SHTNS_INC)
    list(APPEND INCLUDE_DIRS ${_SHTNS_INC})
  else(_SHTNS_INC)
    message(FATAL_ERROR " Could not detect the SHTNS include directory.")
  endif(_SHTNS_INC)
else(_SHTNS_LIB)
  message(FATAL_ERROR " Could not detect the SHTNS library directory.")
endif(_SHTNS_LIB)
unset(_SHTNS_LIB CACHE)

# find FFTW libraries
find_library(_FFTW_LIB fftw3 HINTS ${FFTW_LIB} $ENV{FFTW_LIB})
if(_FFTW_LIB)
  get_filename_component(_FFTW_DIR ${_FFTW_LIB} PATH CACHE)
  list(APPEND LIBRARY_DIRS ${_FFTW_DIR})
  list(APPEND LIBRARIES "${_FFTW_LIB}")

  find_library(_FFTW_LIB2 fftw3_omp HINTS ${_FFTW_DIR})
  if(_FFTW_LIB2)
    list(APPEND LIBRARIES "${_FFTW_LIB2}")
  endif(_FFTW_LIB2)

  # find FFTW include-directory
  find_path(_FFTW_INC fftw3.h HINTS $ENV{FFTW_INC} ${_FFTW_DIR}/../include)
  if(_FFTW_INC)
    list(APPEND INCLUDE_DIRS ${_FFTW_INC})
  else(_FFTW_INC)
    message(FATAL_ERROR " Could not detect the FFTW include directory.")
  endif(_FFTW_INC)
else(_FFTW_LIB)
  message(FATAL_ERROR " Could not detect the FFTW library directory.")
endif(_FFTW_LIB)
unset(_FFTW_LIB CACHE)

find_file(_EIGEN Eigen HINTS $ENV{EIGEN_INC}/Eigen $ENV{EIGEN3_INC}/eigen3/Eigen /usr/include/eigen3/Eigen)
if(_EIGEN)
    get_filename_component(_EIGEN_DIR ${_EIGEN} PATH CACHE)
    list(APPEND INCLUDE_DIRS ${_EIGEN_DIR}/..)
else(_EIGEN)
  message(FATAL_ERROR " Could not detect the Eigen include directory.")
endif(_EIGEN)

message(STATUS "LIBRARIES: ${LIBRARIES}")
message(STATUS "INCLUDE_DIRS: ${INCLUDE_DIRS}")

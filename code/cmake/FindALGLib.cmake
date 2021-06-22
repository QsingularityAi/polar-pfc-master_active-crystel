# find ALGLib libraries
find_library(_ALGLIB_LIB alglib
  HINTS ${ALGLIB_LIB} $ENV{ALGLIB_LIB})

if(_ALGLIB_LIB)
  get_filename_component(_ALGLIB_DIR ${_ALGLIB_LIB} PATH CACHE)
endif()

# find ALGLib include directory
find_path(_ALGLIB_INC alglibmisc.h
  HINTS ${ALGLIB_INC} $ENV{ALGLIB_INC} ${_ALGLIB_DIR}/../../include/alglib
  PATH_SUFFIXES include libalglib include/libalglib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ALGLib _ALGLIB_LIB _ALGLIB_INC)
mark_as_advanced(_ALGLIB_LIB _ALGLIB_INC _ALGLIB_DIR)

if (ALGLib_FOUND AND NOT TARGET ALGLib::ALGLib)
  add_library(ALGLib::ALGLib UNKNOWN IMPORTED)
  set_target_properties(ALGLib::ALGLib PROPERTIES
    IMPORTED_LOCATION ${_ALGLIB_LIB}
    INTERFACE_INCLUDE_DIRECTORIES ${_ALGLIB_INC}
  )
endif()
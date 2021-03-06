# Pastix lib requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.
# Pastix requires SCOTCH or METIS (partitioning and reordering tools) as well

IF(DEFINED ENV{PASTIX_ROOT})
   SET(PASTIX_ROOT $ENV{PASTIX_ROOT} CACHE PATH "pastix location")
ELSE()
   SET(PASTIX_ROOT /usr/local CACHE PATH "pastix location")
ENDIF()

FIND_PATH(PASTIX_INCLUDE_DIRS
	    NAMES pastix_fortran.h
	    HINTS ${PASTIX_ROOT}
	    PATH_SUFFIXES include Include INCLUDE
	    DOC "PATH TO pastix_fortran.h")

FIND_LIBRARY(PASTIX_LIBRARY NAMES pastix
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libpastix.a")

FIND_LIBRARY(PASTIX_MATRIX_DRIVER_LIBRARY NAMES matrix_driver
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libmatrix_driver.a")

FIND_LIBRARY(HWLOC_LIBRARY NAMES hwloc
		 HINTS $ENV{HWLOC_HOME} ${PASTIX_ROOT} /usr/local /opt/local
		 PATH_SUFFIXES lib Lib LIB lib64
		 DOC "PATH TO hwloc library")

SET(PASTIX_LIBRARIES ${PASTIX_LIBRARY};${PASTIX_MATRIX_DRIVER_LIBRARY};${HWLOC_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PASTIX DEFAULT_MSG PASTIX_INCLUDE_DIRS PASTIX_LIBRARIES)

#IF(PASTIX_FOUND)
   MESSAGE(STATUS "PASTIX_INCLUDE_DIRS:${PASTIX_INCLUDE_DIRS}")
   MESSAGE(STATUS "PASTIX_LIBRARIES:${PASTIX_LIBRARIES}")
#ENDIF(PASTIX_FOUND)

FIND_PATH(MURGE_INCLUDE_DIRS
	    NAMES murge.inc
	    HINTS ${PASTIX_ROOT}
	    PATH_SUFFIXES include Include INCLUDE
	    DOC "PATH TO murge.inc")

FIND_LIBRARY(MURGE_LIBRARIES
		 NAMES pastix_murge 
		 HINTS ${PASTIX_ROOT}
		 PATH_SUFFIXES lib Lib LIB
		 DOC "PATH TO libpastix_murge.a")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MURGE DEFAULT_MSG MURGE_INCLUDE_DIRS MURGE_LIBRARIES)

MARK_AS_ADVANCED(PASTIX_INCLUDE_DIRS
                 PASTIX_LIBRARY
                 PASTIX_MATRIX_DRIVER_LIBRARY
                 PASTIX_LIBRARIES
                 HWLOC_LIBRARY
                 MURGE_INCLUDE_DIRS
                 MURGE_LIBRARIES)
                 

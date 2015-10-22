# - Try to find Jorek 
# Once done this will define
#
#  JOREK_DIR          - Install Directory for Jorek 
#  JOREK_FOUND        - system has Jorek
#  JOREK_INCLUDES     - Jorek include directories
#  JOREK_LIBRARIES    - Link these to use Jorek
#  The following variables are TODO
#  JOREK_COMPILER     - Compiler used by Jorek, helpful to find a compatible MPI
#  JOREK_DEFINITIONS  - Compiler switches for using Jorek
#  JOREK_MPIEXEC      - Executable for running MPI programs
#  JOREK_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
#
#  Usage:
#  find_package(JOREK)                  
#
# Setting these changes the behavior of the search
#  JOREK_DIR - directory in which PETSc resides
#  The following variable is TODO
#  JOREK_ARCH - build architecture
#

SET(JOREK_DIR "$ENV{JOREK_DIR}")
SET(JOREK_INCLUDES "${JOREK_DIR}/include")
SET(JOREK_LIB_DIRS "${JOREK_DIR}/lib")

FIND_LIBRARY(JOREK_LIB_TRACELOG NAMES tracelog
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_GLOB NAMES jorek_glob
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_PARAM NAMES jorek_param
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_CONTEXT NAMES jorek_context
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_SPM NAMES spm
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_NUMBERING NAMES numbering
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_MESH NAMES mesh
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_COORDINATES NAMES coordinates
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_BASIS NAMES basis
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_ASSEMBLY NAMES assembly
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_EVALUATOR NAMES evaluator
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_OUTPUT NAMES output
	                  HINTS ${JOREK_LIB_DIRS})

FIND_LIBRARY(JOREK_LIB_BC NAMES boundary_conditions
	                  HINTS ${JOREK_LIB_DIRS})

SET(JOREK_LIBRARIES 
	${JOREK_LIB_TRACELOG}
	${JOREK_LIB_GLOB}
	${JOREK_LIB_PARAM}
	${JOREK_LIB_CONTEXT}
	${JOREK_LIB_SPM}
	${JOREK_LIB_NUMBERING}
	${JOREK_LIB_MESH}
	${JOREK_LIB_COORDINATES}
	${JOREK_LIB_BASIS}
	${JOREK_LIB_ASSEMBLY}
	${JOREK_LIB_EVALUATOR}
	${JOREK_LIB_OUTPUT}
	${JOREK_LIB_BC}
	)

SET (JOREK_INCLUDES  ${JOREK_INCLUDES} CACHE STRING "Jorek include path" FORCE)
SET (JOREK_LIBRARIES ${JOREK_LIBRARIES} CACHE STRING "Jorek libraries" FORCE)
MARK_AS_ADVANCED(JOREK_INCLUDES JOREK_LIBRARIES)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args (JOREK
	"Jorek could not be found.  Be sure to set JOREK_DIR."
	JOREK_INCLUDES JOREK_LIBRARIES)

cmake_minimum_required (VERSION 2.6)
project (Seqan_Tests)

# Some tests need OpenMP.
find_package (OpenMP QUIET)

################################################################################
# For the tests, we switch SeqAn into testing mode.
################################################################################

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${SAMTOOLS_CXX_FLAGS} -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=1")
# Also enable flags for coverage analysis.
if (MODEL STREQUAL "NightlyCoverage")
  if (CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
  endif (CMAKE_COMPILER_IS_GNUCXX)
endif (MODEL STREQUAL "NightlyCoverage")
if (MODEL STREQUAL "ExperimentalCoverage")
  if (CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
  endif (CMAKE_COMPILER_IS_GNUCXX)
endif (MODEL STREQUAL "ExperimentalCoverage")

################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)

################################################################################
# Define which tests are to link against libbam.
################################################################################

MACRO(LIST_CONTAINS var value)
  SET(${var})
  FOREACH (value2 ${ARGN})
    IF (${value} STREQUAL ${value2})
      SET(${var} TRUE)
    ENDIF (${value} STREQUAL ${value2})
  ENDFOREACH (value2)
ENDMACRO(LIST_CONTAINS)

set(SEQAN_LINK_LIBBAM ${SEQAN_LINK_LIBBAM} test_store)

################################################################################
# Tests
################################################################################

add_custom_target(tests)

#-----------------------------------------------------------------------------
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.
INCLUDE(CTest)
ENABLE_TESTING()
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

get_filename_component (SEQAN_TESTS_ABS ${SEQAN_LIBRARY}/../tests ABSOLUTE)
#if (EXISTS SEQAN_TESTS_ABS)
	file(GLOB SEQAN_TESTS ${SEQAN_TESTS_ABS}/[A-z]*)
	foreach (TEST_DIR ${SEQAN_TESTS})
		file (RELATIVE_PATH TEST_DIR_REL ${SEQAN_TESTS_ABS} ${TEST_DIR})
		file (GLOB TEST_SOURCES ${TEST_DIR}/[A-z]*.h ${TEST_DIR}/test_[A-z]*.cpp)
		list (SORT TEST_SOURCES)
		if (TEST_SOURCES)
			add_executable (test_${TEST_DIR_REL} ${TEST_SOURCES})
			add_dependencies (test_${TEST_DIR_REL} ${SEQAN_TARGET})
	    add_dependencies (tests test_${TEST_DIR_REL})

    # Conditionally link against libbam.
    if (SAMTOOLS_FOUND)
      LIST_CONTAINS(contains test_${TEST_DIR_REL} ${SEQAN_LINK_LIBBAM})
      if (contains)
        if (ZLIB_FOUND)
          target_link_libraries(test_${TEST_DIR_REL} ${SAMTOOLS_LIBRARIES} ${ZLIB_LIBRARIES})
        else (ZLIB_FOUND)
          target_link_libraries(test_${TEST_DIR_REL} ${SAMTOOLS_LIBRARIES})
        endif (ZLIB_FOUND)
      endif (contains)
    endif (SAMTOOLS_FOUND)

      # Create a variable with the paths to the tests for memory checks.
      # TODO(holtgrew): For MSVC the paths are different, but does not matter as not run on Windows yet.
      if (MODEL MATCHES ".*Mem.*Check")
        set(TEST_PATH ${SEQAN_BINARY_ROOT}/projects/library/cmake/tests/)
      endif (MODEL MATCHES ".*Mem.*Check")
			add_test (test_${TEST_DIR_REL} ${TEST_PATH}test_${TEST_DIR_REL})
			if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
				target_link_libraries (test_${TEST_DIR_REL} rt)
			endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		endif (TEST_SOURCES)
	endforeach (TEST_DIR)
#endif (EXISTS SEQAN_TESTS_ABS)

set (SEQAN_TESTS_APPS_DIR "${SEQAN_TESTS_ABS}/apps")

file(GLOB SEQAN_TESTS ${SEQAN_TESTS_APPS_DIR}/[A-z]*)
foreach (TEST_DIR ${SEQAN_TESTS})
  file (RELATIVE_PATH TEST_DIR_REL ${SEQAN_TESTS_APPS_DIR} ${TEST_DIR})
  if (MODEL MATCHES ".*Mem.*Check")
    add_test(test_apps_${TEST_DIR_REL} ${PYTHON_EXECUTABLE} ${TEST_DIR}/run_tests.py --valgrind ${SEQAN_SOURCE_ROOT} ${SEQAN_BINARY_ROOT})
  else (MODEL MATCHES ".*Mem.*Check")
    add_test(test_apps_${TEST_DIR_REL} ${PYTHON_EXECUTABLE} ${TEST_DIR}/run_tests.py ${SEQAN_SOURCE_ROOT} ${SEQAN_BINARY_ROOT})
  endif (MODEL MATCHES ".*Mem.*Check")
endforeach (TEST_DIR)


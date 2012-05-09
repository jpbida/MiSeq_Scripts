cmake_minimum_required (VERSION 2.6)
project (Seqan_Apps)

################################################################################
# For the apps, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
################################################################################

# Use CMake's package location functionality to find the directories of
find_package (OpenMP QUIET)
find_package (TR1)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${SAMTOOLS_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=0")
set (CMAKE_CXX_FLAGS_RELDEBUG "${CMAKE_CXX_FLAGS_RELDEBUG} -DSEQAN_ENABLE_DEBUG=0")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

################################################################################
# Set path variables
################################################################################

# Make sure the compiler can find include files from the SeqAn.
include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)

################################################################################
# Applications
################################################################################

# TODO(holtgrew): Disable apps if dependencies are missing.  This would mean removing the automatic generation of targets but this invites errors anyway.

# pseudo-target to build all apps
add_custom_target (apps)

# Find targets automatically
get_filename_component (SEQAN_APPS_ABS ${SEQAN_LIBRARY}/apps ABSOLUTE)
file (GLOB SEQAN_APPS ${SEQAN_APPS_ABS}/[A-z]*)

file (GLOB RAZERS_HEADERS ${SEQAN_LIBRARY}/apps/razers/[A-z]*.h)
file (GLOB RAZERS2_HEADERS ${SEQAN_LIBRARY}/apps/razers2/[A-z]*.h)
file (GLOB SEQAN_TCOFFEE_HEADERS ${SEQAN_LIBRARY}/apps/seqan_tcoffee/[A-z]*.h)

# Manually set predefined sets of source files for some targets
set (APP_razers_SOURCES ${SEQAN_LIBRARY}/apps/razers/razers.cpp ${RAZERS_HEADERS})
set (APP_razers2_SOURCES ${SEQAN_LIBRARY}/apps/razers2/razers.cpp ${RAZERS2_HEADERS})
set (APP_micro_razers_SOURCES ${SEQAN_LIBRARY}/apps/micro_razers/micro_razers.cpp ${RAZERS_HEADERS})
set (APP_param_chooser_SOURCES ${SEQAN_LIBRARY}/apps/param_chooser/param_chooser.cpp ${RAZERS_HEADERS})

################################################################################
# Auto-generate headers.
################################################################################

foreach (APP_DIR ${SEQAN_APPS})
  if (EXISTS ${APP_DIR}/config.h.in)
    configure_file (${APP_DIR}/config.h.in
                    ${APP_DIR}/config.h)
  endif (EXISTS ${APP_DIR}/config.h.in)
endforeach (APP_DIR ${SEQAN_APPS})

################################################################################
# Build applications.
################################################################################

# Make sure the compiler can find object files from the TBB library.
link_directories (${SEQAN_LIBRARY}/lib)

foreach (APP_DIR ${SEQAN_APPS})
	file (RELATIVE_PATH APP_DIR_REL ${SEQAN_APPS_ABS} ${APP_DIR})

	# if there is no predefined set of files for this target
	# find the source files automatically
	if (NOT APP_${APP_DIR_REL}_SOURCES)
		file (GLOB APP_${APP_DIR_REL}_SOURCES ${APP_DIR}/[A-z]*.cpp ${APP_DIR}/[A-z]*.h ${APP_DIR}/[A-z]*.cu)
	endif (NOT APP_${APP_DIR_REL}_SOURCES)

	foreach (SRC ${APP_${APP_DIR_REL}_SOURCES})
		if (NOT EXISTS ${SRC})
#			message ("WARNING ${SRC} not found. Disabling ${APP_DIR_REL}.")
			set (APP_${APP_DIR_REL}_NOBUILD TRUE)
		endif (NOT EXISTS ${SRC})
	endforeach (SRC)

	# add target, dependencies and linked libraries
	set (SAVED_FLAGS ${CMAKE_CXX_FLAGS})
	if (APP_${APP_DIR_REL}_SOURCES)
		if (APP_${APP_DIR_REL}_NOBUILD)
			message ("Not building ${APP_DIR_REL} -- it is deactivated.")
		else (APP_${APP_DIR_REL}_NOBUILD)
			list (SORT APP_${APP_DIR_REL}_SOURCES)
			  add_executable (${APP_DIR_REL} ${APP_${APP_DIR_REL}_SOURCES})
			  install(TARGETS ${APP_DIR_REL}
			          DESTINATION bin
			          COMPONENT apps)
			add_dependencies (${APP_DIR_REL} ${SEQAN_TARGET})
			add_dependencies (apps ${APP_DIR_REL})
            # Link against libz everywhere, used for reading BAM, for example.
            if (ZLIB_FOUND)
			  target_link_libraries (${APP_DIR_REL} z)
            endif (ZLIB_FOUND)
            # Link against librt on Linux for real-time clock.
			if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
				target_link_libraries (${APP_DIR_REL} rt)
			endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
			set (CMAKE_CXX_FLAGS ${SAVED_FLAGS})
		endif (APP_${APP_DIR_REL}_NOBUILD)
	endif (APP_${APP_DIR_REL}_SOURCES)
endforeach (APP_DIR)

# Install RazerS gapped parameters.
install(DIRECTORY ${SEQAN_LIBRARY}/apps/razers2/gapped_params
        DESTINATION bin)

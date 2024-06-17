cmake_minimum_required(VERSION 3.1)
project(MercuryDPM C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)

# Trunk will become 0.14.alpha
set(MercuryDPM_VERSION_MAJOR 1)
set(MercuryDPM_VERSION_MINOR 0)

set(REVISION, "0")
set(REPOSITORY_URL "")

find_package(Git)
if (EXISTS ${PROJECT_SOURCE_DIR}/.git)
	if (Git_FOUND)
		execute_process(
				COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
				WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
			OUTPUT_VARIABLE REVISION
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		# For Git it doesn't make sense to specify a URL since the
		# commit might not have been pushed to a server. (Using
		# `git remote -v` doesn't help if the repository has more than
		# one remote set up.)
		set(REPOSITORY_URL "https://bitbucket.org/mercurydpm/mercurydpm.git")
		#message("Git: commit is ${REVISION}")
	else()
		message("The project is in a Git repository, but the git executable is not available")
	endif()
else()
	message("This MercuryDPM version is not in Git")
	find_package(Subversion)
	if(EXISTS ${PROJECT_SOURCE_DIR}/.svn)
		if(Subversion_FOUND)
			Subversion_WC_INFO(${PROJECT_SOURCE_DIR} SVN)
			set(REPOSITORY_URL ${SVN_WC_URL})
			set(REVISION ${Project_WC_REVISION})
		endif()
	else()
		message("This MercuryDPM version in not in svn")
	endif()
endif()

set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

enable_testing()

if (UNIX)
    add_definitions( -DUNIX )
endif(UNIX)

if (WIN32)
    message("Windows operating system")
    add_definitions( -DWINDOWS )
endif(WIN32)

if (CYGWIN)
    message("You are using CYGWIN; which is no longer supported. However the code should still work")
    add_definitions(-DCYGWIN)
endif (CYGWIN)

#If user has not given build type set it to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE
            STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#If cMake builds in debug mode add the MercuryDPM_DEBUG flag
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DMERCURYDPM_DEBUG)
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

include(CMakeModules/MercuryDPM_MPI.cmake)

include(CMakeModules/MercuryDPM_OMP.cmake)

option(MercuryDPM_TRIANGLE_WALL_CORRECTION "Smooth transition between triangle walls" ON)
message(STATUS "MercuryDPM_TRIANGLE_WALL_CORRECTION" ${MercuryDPM_TRIANGLE_WALL_CORRECTION})
message(STATUS "MercuryDPM_USE_OpenMP" ${MercuryDPM_USE_OpenMP})
if(MercuryDPM_TRIANGLE_WALL_CORRECTION)
	if(MercuryDPM_USE_OpenMP)
		message(STATUS "You are using MercuryDPM_USE_OpenMP; this is incompatible with MercuryDPM_TRIANGLE_WALL_CORRECTION. Therefore MercuryDPM_TRIANGLE_WALL_CORRECTION will be disabled")
	else ()
		add_definitions( -DMERCURYDPM_TRIANGLE_WALL_CORRECTION )
	endif ()
endif()

include(CMakeModules/MercuryDPM_Logger.cmake)

include(CMakeModules/MercuryDPM_Backtrace.cmake)

include(CMakeModules/MercuryDPM_CppFeatureCheck.cmake)

include(CMakeModules/MercuryDPM_Documentation.cmake)

include(CMakeModules/MercuryDPM_Xballs.cmake)

#This section generates a configuration file that contains the CMake variables#
#############################################
configure_file(Configuration/CMakeDefinitions.cc
		Kernel/CMakeDefinitions.cc)
configure_file(Configuration/CMakeDefinitions.h
		Kernel/CMakeDefinitions.h)

option(MercuryDPM_KERNEL_ONLY "Only the kernel files are built, no drivers, no tests." OFF)
mark_as_advanced(MercuryDPM_KERNEL_ONLY)

add_subdirectory(Kernel)

if (NOT MercuryDPM_KERNEL_ONLY)
	add_subdirectory(Tools)
	add_subdirectory(XBalls)
endif()

# Determine whether to couple MercuryDPM with oomph-lib or not
option(MercuryDPM_OOMPH_COUPLING "Couple MercuryDPM with oomph-lib" OFF)
# Couple with oomph-lib if OOMPH_COUPLING is on
include(CMakeModules/OomphCoupling.cmake)

# Determine whether to couple MercuryDPM with Eigen or not
option(MercuryDPM_EIGEN "Couple MercuryDPM with Eigen" OFF)
# Couple with Eigen if EIGEN is on
include(CMakeModules/MercuryDPM_Eigen.cmake)

if (NOT MercuryDPM_KERNEL_ONLY)
	# this is done after oomph and Eigen is included since it needs the includes defined there
	add_subdirectory(Drivers)

	#Add the command fullTest, which rebuilds the cache, cleans the project and then compiles and tests everything.
	add_custom_target(fullTest
		COMMAND make rebuild_cache
		COMMAND make clean
		COMMAND make -j 4
		COMMAND make test -j 4
		COMMENT "This will clean the project, rebuild and test." )
endif()

#Add the command update, which rebuilds the cache and compiles the codes (no test).
add_custom_target(update
	COMMAND make rebuild_cache
	COMMAND make -j 4
	COMMENT "Quick rebuild" )

#This section is about the build of documentation#
##################################################

option(MercuryDPM_BUILD_USER_DIR "Compile USER DIR in testing" ON)
if(MercuryDPM_KERNEL_ONLY AND MercuryDPM_BUILD_USER_DIR)
	message(WARNING "MercuryDPM_KERNEL_ONLY overwrites MercuryDPM_BUILD_USER_DIR")
	set(MercuryDPM_BUILD_USER_DIR OFF)
endif()

if (NOT MercuryDPM_KERNEL_ONLY)
	###### These needs to be fixed as the new scripts line above will always create this, so the Python problem is back
	################################################
	find_package(PythonInterp 3 REQUIRED)
	find_package(PythonLibs 3)

	# this line causes problems when you have multiple builds
	FILE(COPY Configuration/fpdiff.py DESTINATION ${PROJECT_SOURCE_DIR}/Scripts/)

	#find_program(PYTHON_EXEC "python"
	#	DOC "Location of python executable to use")
endif()

#Do nothing if coupling is turned off.
if(NOT MercuryDPM_OOMPH_COUPLING)
    #message(STATUS "Coupling with oomph-lib is disabled")
    return()
endif()

#check if git is installed
if (NOT Git_FOUND)
    message(FATAL_ERROR "The option you have chosen requires git and you do not have this installed. Please install")
endif()

# Clone oomph-lib if has not been cloned before
set(OOMPH_DIR ${PROJECT_SOURCE_DIR}/oomph-lib)
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND git submodule init ${OOMPH_DIR})
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND git submodule update  ${OOMPH_DIR})
if(NOT EXISTS ${OOMPH_DIR}/src)
    message(FATAL_ERROR "git clone failed. If this problem persists you can manually clone oomph-lib by running: \n   git submodule init\n   git submodule update")
endif()

#add CMakeCache option
option(OOMPH_HAS_MPI "Whether MPI is installed" OFF)
if (OOMPH_HAS_MPI)
    add_definitions( -DOOMPH_HAS_MPI )
    FIND_PACKAGE(MPI REQUIRED)
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "The option you have chosen requires mpi and you do not have this installed. Please install")
    endif()
    message(STATUS "MPI found: ${MPI_CXX_INCLUDE_PATH}")
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_definitions( -DOOMPH_HAS_MPI )
    # add libraries that are only needed in MPI mode
    set(MPI_CXX_LIBRARIES oomph_superlu_dist_3.0 oomph_metis_from_parmetis_3.1.1 ${MPI_CXX_LIBRARIES})
    add_definitions( -DUSING_OOMPH_SUPERLU_DIST )
else()
    set(MPI_CXX_LIBRARIES "")
endif()

#add CMakeCache options
option(MUMPS_INCLUDE_DIR "Whether MUMPS is installed" "")
if (MPI_FOUND AND MUMPS_INCLUDE_DIR)
    FIND_PACKAGE(MUMPS)
    if (MUMPS_FOUND)
        add_definitions( -DOOMPH_HAS_MUMPS )
        message(STATUS "MUMPS found: ${MUMPS_INCLUDE_DIR}")
        #add_library(MUMPS_LIBRARY)
        include_directories(${MUMPS_INCLUDE_DIR})
    else()
        message(ERROR "MUMPS not found; please set OOMPH_HAS_MUMPS=OFF")
    endif()
endif()

# Choose whether to compile with cmake or use prebuilt oomph
option(OOMPH_CMAKE "Use cmake to couple oomph-lib (if coupling turned on)" ON)
if (OOMPH_CMAKE)
    # copy CMakeFiles, add oomph directory to compilation
    message(STATUS "Compiling oomph with CMake")
    FILE(COPY CMakeModules/oomph/CMakeLists.txt DESTINATION ${OOMPH_DIR})
    FILE(COPY CMakeModules/oomph/src/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src)
    FILE(COPY CMakeModules/oomph/src/generic/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src/generic)
    FILE(COPY CMakeModules/oomph/src/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src/solid)
    FILE(COPY CMakeModules/oomph/external_src/CMakeLists.txt DESTINATION ${OOMPH_DIR}/external_src)
    FILE(COPY CMakeModules/oomph/demo_drivers/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers)
    FILE(COPY CMakeModules/oomph/demo_drivers/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers/solid)
    FILE(COPY CMakeModules/oomph/demo_drivers/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers/solid)
    if (${APPLE})
        message (STATUS "MAC OS X: making dummy fpu_control.h")
        IF(NOT EXISTS ${OOMPH_DIR}/external_src/oomph_triangle/fpu_control.h)
            FILE(COPY ${OOMPH_DIR}/external_src/oomph_triangle/dummy_fpu_control.h DESTINATION ${OOMPH_DIR})
            FILE(RENAME ${OOMPH_DIR}/dummy_fpu_control.h ${OOMPH_DIR}/external_src/oomph_triangle/fpu_control.h)
        ENDIF()
        # silence some oomph-warnings
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration -Wformat-extra-args -Wno-parentheses -Wno-implicit-int -Wno-format-security ")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template -Wno-instantiation-after-specialization")
    endif ()
    add_subdirectory(${OOMPH_DIR})
else()
    # Otherwise, we assume it is already compiled
    message(STATUS "Coupling to external oomph directory" ${OOMPH_DIR})
    #check if the oomph-lib Src directory exists
    if (NOT IS_DIRECTORY ${OOMPH_DIR}/src)
        message(FATAL_ERROR "${OOMPH_DIR}/src does not exist.\n $ENV{HOME} Set OOMPH_DIR to the directory where the src/ and build/lib/ folders of oomphlib reside")
    endif()
    #check if the oomph-lib standard libraries have been built
    if (NOT IS_DIRECTORY ${OOMPH_DIR}/build/lib)
        message(FATAL_ERROR "${OOMPH_DIR}/build/lib does not exist.\n Please compile oomph-lib first")
    endif()
    # include the folder where the oomph-lib libraries reside
    # (so we can link against the libraries, see e.g. Drivers/OomphlibCoupling/CMakeLists.txt)
    link_directories(${OOMPH_DIR}/build/lib)
endif()

# include the folders where the oomph-lib Src files reside
# (so you can write e.g. #include "mesh.h")
include_directories(${OOMPH_DIR}/src ${OOMPH_DIR}/src/poisson ${OOMPH_DIR}/src/generic ${OOMPH_DIR}/src/solid ${OOMPH_DIR}/src/constitutive ${OOMPH_DIR}/external_src)


# link some essential external libraries to generic as a minimal oomph library (adding daxpy.f just serves as a dummy source file)
add_library(oomphBase STATIC ${CMAKE_SOURCE_DIR}/Kernel/Math/daxpy.f)
# note libraries need to be ordered such that the left-most can depend on the right most, but not the other way round, i.e. simplest libraries rightmost/last.
target_link_libraries(oomphBase generic ${MPI_CXX_LIBRARIES} ${MUMPS_LIBRARIES} oomph_superlu_4.3 oomph_flapack oomph_arpack oomph_blas oomph_lapack)

# build a smaller library for solid problems
add_library(oomphSolid STATIC ${CMAKE_SOURCE_DIR}/Kernel/Math/daxpy.f)
target_link_libraries(oomphSolid constitutive meshes solid oomphBase)

# the full oomph library
add_library(oomph STATIC ${CMAKE_SOURCE_DIR}/Kernel/Math/daxpy.f)
target_link_libraries(oomph steady_axisym_advection_diffusion young_laplace advection_diffusion advection_diffusion_reaction axisym_advection_diffusion axisym_foeppl_von_karman axisym_linear_elasticity axisym_navier_stokes axisym_poroelasticity axisym_spherical_solid beam biharmonic constitutive darcy fluid_interface flux_transport foeppl_von_karman fourier_decomposed_helmholtz generalised_newtonian_axisym_navier_stokes generalised_newtonian_navier_stokes helmholtz linear_elasticity linear_wave linearised_navier_stokes linearised_axisym_navier_stokes mesh_smoothing meshes multi_physics navier_stokes ode poisson polar_navier_stokes poroelasticity rigid_body shell solid spherical_advection_diffusion spherical_navier_stokes  time_harmonic_fourier_decomposed_linear_elasticity time_harmonic_linear_elasticity unsteady_heat womersley oomph_hsl oomph_crbond_bessel oomph_triangle oomph_tetgen oomphBase)
# missing:  reynolds_averaged_navier_stokes

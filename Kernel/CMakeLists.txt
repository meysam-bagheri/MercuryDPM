#Include directories where the source files and the libraries are stored
include_directories(
        ${MercuryDPM_SOURCE_DIR}/Kernel
        ${MercuryDPM_BINARY_DIR}/Kernel
)
if (MercuryDPM_USE_MPI)
        include_directories(${MPI_CXX_INCLUDE_PATH})
endif ()

#For the Fortran-based linear algebra
add_subdirectory(Math)

#Compile all .cc in the directories below, and add them to the library
file(GLOB MercuryDPM_SRC
        BaseObject.cc
        BaseInteractable.cc
        BaseCluster.cc
        Logger.cc
        File.cc
        DPMBase.cc
        MpiContainer.cc
        MpiDataClass.cc
        HGrid.cc
        MercuryBase.cc
        Mercury2D.cc
        Mercury3D.cc
        HGridOptimiser.cc
        BoundaryHandler.cc
        PeriodicBoundaryHandler.cc
        InteractionHandler.cc
        ParticleHandler.cc
        SpeciesHandler.cc
        WallHandler.cc
        DomainHandler.cc
        Domain.cc
        BinaryReader.cc
        ${PROJECT_BINARY_DIR}/Kernel/CMakeDefinitions.cc

        Boundaries/*.cc
        CG/*.cc
        CG/*/*.cc
        FMM/*.cc
        Math/*.cc
        Math/*/*.cc
        Mixins/*.cc
        Particles/*.cc
        Strings/*.cc
        VTKWriter/*.cc
        Walls/*.cc
        Walls/*/*.cc
        Interactions/*.cc
        Interactions/*/*.cc
        Species/*.cc
        Species/*/*.cc
        Nurbs/*.cc
        Tools/*.cc
        Controllers/*.cc
        Helpers/*.cc
        )

#Basic library
add_library(MercuryBase STATIC ${MercuryDPM_SRC})
target_link_libraries(MercuryBase LinearAlgebra)

if (MercuryDPM_BACKTRACE_DEMANGLE)
        if (MercuryDPM_BACKTRACE_LINK_DL)
                target_link_libraries(MercuryBase dl)
        endif ()
        if (MercuryDPM_BACKTRACE_LINK_EXECINFO)
                target_link_libraries(MercuryBase execinfo)
        endif ()
endif ()

if (MercuryDPM_USE_MPI)
        set(DEP_COMPILE_FLAGS "${DEP_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
        set(DEP_LINKER_FLAGS "${DEP_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
        target_link_libraries(MercuryBase "${MPI_LIBRARIES}")
endif ()

#Required for Chute drivers
add_library(Chute STATIC
        ChuteBottom.cc
        ChuteWithHopper.cc
        Chute.cc
        CurvyChute.cc
        )

target_link_libraries(Chute MercuryBase)

#This file is included in all CMakeLists.txt files in the Driver directory

#Define includes for compiling the Driver codes
#The Mercury_BINARY_DIR is added here to enable multiple builds in a single file system
#The CMakeDefinitions.cc and CMakeDefinitions.h file are copied into the binary folder
include_directories(${MercuryDPM_SOURCE_DIR}/Kernel
		${MercuryDPM_BINARY_DIR}/Kernel)

#Part 2 : Make run test for each of the demo files
##################################################

#For MPI*Test drivers:
#Define a function to extract the number of required processors
function(get_number_of_cores EXECNAME NUMCORES)
	string(FIND ${EXECNAME} "MPI" POS1)
  string(FIND ${EXECNAME} "Test" POS2)
  math(EXPR START "${POS1} + 3")
  math(EXPR LENGTH "${POS2} - ${START}") 
  if (${LENGTH} STREQUAL 0)
    message(FATAL_ERROR "No number of cores specified for ${EXECNAME}. Format is *MPI<number_of_cores>Test.cpp")
  endif()
  string(SUBSTRING ${EXECNAME} ${START} ${LENGTH} NUMCORES )
  set(NUMCORES ${NUMCORES} PARENT_SCOPE)
endfunction()

file(GLOB SELFTESTS "*SelfTest.cpp")
file(GLOB UNITTESTS "*UnitTest.cpp")
file(GLOB MPITESTS  "*MPI*Test.cpp")
#For each demo add a test with the same name

foreach (TEST ${UNITTESTS} ${SELFTESTS})
	get_filename_component(EXECNAME ${TEST} NAME_WE)
	get_filename_component(EXECPATH ${TEST} DIRECTORY)
	file(RELATIVE_PATH EXECDIR ${PROJECT_SOURCE_DIR}/Drivers ${EXECPATH}) 
	string(CONCAT TESTNAME ${EXECDIR} ":" ${EXECNAME})
	#message(STATUS ${TESTNAME})	
	add_test(${TESTNAME} ${EXECNAME})
    #Set a 30 second timeout for each test.
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 80)
endforeach()

if (MercuryDPM_USE_MPI)
	foreach (TEST ${MPITESTS})
        	get_filename_component(EXECNAME ${TEST} NAME_WE)
            get_number_of_cores(${EXECNAME} NUMCORES)
        	add_test(${EXECNAME} mpiexec -n ${NUMCORES} ./${EXECNAME})
            #Set a 30 second timeout for each test.
            set_tests_properties(${EXECNAME} PROPERTIES TIMEOUT 80)
	endforeach()
endif()

#Part 3 : Make tests for each of the SelfTestData files
########################################################


file(GLOB SELFTESTDATAFILES "${CMAKE_CURRENT_SOURCE_DIR}/SelfTestData/*.*")
file(GLOB MPITESTDATAFILES "${CMAKE_CURRENT_SOURCE_DIR}/MPITestData/*.*")

#for each file in the selftest_data folder create a test. Which checks the data against this old data. The actually testing is done my the script self_test.
foreach(TESTFILE ${SELFTESTDATAFILES})
        	get_filename_component(EXECNAME ${TESTFILE} NAME)
		get_filename_component(EXECPATH ${TESTFILE} DIRECTORY)
		file(RELATIVE_PATH EXECDIR ${PROJECT_SOURCE_DIR}/Drivers ${EXECPATH}) 
		string(CONCAT TESTNAME ${EXECDIR} ":" ${EXECNAME})
		#message(STATUS ${TESTNAME})
		add_test(${TESTNAME} ${CMAKE_SOURCE_DIR}/Scripts/self_test ${TESTFILE} ${EXECNAME})
        #Add the newly created files to the clean target
        set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${TESTNAME}")
endforeach()

if (MercuryDPM_USE_MPI)
	foreach(TESTFILE ${MPITESTDATAFILES})
		get_filename_component(TESTNAME ${TESTFILE} NAME)
		get_filename_component(EXECPATH ${TESTFILE} DIRECTORY)
		file(RELATIVE_PATH EXECDIR ${PROJECT_SOURCE_DIR}/Drivers ${EXECPATH}) 
		string(CONCAT TESTNAME ${EXECDIR} ":" ${TESTNAME})
		#message(STATUS ${TESTNAME})	
		add_test(${TESTNAME} ${CMAKE_SOURCE_DIR}/Scripts/self_test ${TESTFILE} ${TESTNAME})
		#Add the newly created files to the clean target
		set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${TESTNAME}")
	endforeach()
endif()



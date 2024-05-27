include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

###########################
#Part 1 : All cpp files will be made to exec files
#########################################
#Collect all the names of the cpp, note at the moment limited to demos, but will be fixed later
file(GLOB CPPFILES  "*.cpp")
#for every cpp found
foreach(CPPFILE ${CPPFILES})
	#extract the actually file name
	get_filename_component(FILENAME ${CPPFILE} NAME)
	#extract the filename minus the cpp. This will be the name of exe file
	get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
	#Make the exe
	add_executable(${EXECNAME} ${FILENAME})
    #All cpp folder and linked against MercuryBase
	target_link_libraries(${EXECNAME} MercuryBase)

	#Setup scripts directory
	########################################################

	#First create the directory if does not exist this stops the dir does not exist error in the next line
	add_custom_command(TARGET ${EXECNAME} PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory  ${CMAKE_CURRENT_BINARY_DIR}/Scripts/)


	#Copy and copy back
	add_custom_command(TARGET ${EXECNAME} POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_BINARY_DIR}/Scripts/ ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/
	COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/ ${CMAKE_CURRENT_BINARY_DIR}/Scripts/)

endforeach()

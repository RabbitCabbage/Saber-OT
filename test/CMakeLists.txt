# add_compile_options(-fsanitize=address)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# add_link_options(-fsanitize=address)

option(ENABLE_SABER "Enable SABER" OFF)
option(ENABLE_KYBER "Enable KYBER" OFF)

#Testing macro
macro (add_test_executable_with_lib _name libs)
	add_executable(test_${_name} "${_name}.cpp")
    
    if(ENABLE_SABER)
		target_link_libraries(test_${_name} ${EMP-TOOL_LIBRARIES})
        target_compile_definitions(test_${_name} PRIVATE SABER)
        target_include_directories(test_${_name} PRIVATE ${CMAKE_SOURCE_DIR}/emp-ot/saber)
    endif()
    
    if(ENABLE_KYBER)
		ADD_SUBDIRECTORY(KyberOT)
		target_link_libraries(test_${_name} ${EMP-TOOL_LIBRARIES} KyberOT)
        target_compile_definitions(test_${_name} PRIVATE KYBER)
        target_include_directories(test_${_name} PRIVATE ${CMAKE_SOURCE_DIR}/emp-ot/kyber)
    endif()

endmacro()

macro (add_test_case _name)
	add_test_executable_with_lib(${_name} "")
  	add_test(NAME ${_name} COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_run _name)
	add_test_executable_with_lib(${_name} "")
	add_test(NAME ${_name} COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

endmacro()

# Test cases
# add_test_case_with_run(ot)
# add_test_case_with_run(ferret)
# add_test_case(bench_lpn)
if(ENABLE_SABER)
    add_test_case(saber)
endif()

if(ENABLE_KYBER)
    add_test_case(mr_kyber)
endif()
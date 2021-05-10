set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(perf_rlbwt perf_rlbwt.cpp)

add_executable(perf_rlbwt_with_monitor perf_rlbwt.cpp)
set_target_properties(perf_rlbwt_with_monitor PROPERTIES COMPILE_DEFINITIONS "ENABLE_STAT_MONITOR")

set(TEST_OPTIONS
    "LFIG_NAIVE"
    "GLFIG_NAIVE_4"
    "GLFIG_NAIVE_16"
    "GLFIG_SERIALIZED_4"
    "GLFIG_SERIALIZED_16"
    )
foreach(TEST_OPTION ${TEST_OPTIONS})
    set(TEST_SRC_NAME test_${TEST_OPTION})
    add_executable(${TEST_SRC_NAME} test_rlbwt.cpp)
    set_target_properties(${TEST_SRC_NAME} PROPERTIES COMPILE_DEFINITIONS ${TEST_OPTION})
    add_test(NAME ${TEST_SRC_NAME} COMMAND ${TEST_SRC_NAME} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endforeach(TEST_OPTION)

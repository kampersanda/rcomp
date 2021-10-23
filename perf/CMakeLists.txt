set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/perf)

add_executable(perf_rlbwt perf_rlbwt.cpp)
add_executable(perf_rindex perf_rindex.cpp)

# add_executable(perf_rlbwt_with_monitor perf_rlbwt.cpp)
# set_target_properties(perf_rlbwt_with_monitor PROPERTIES COMPILE_DEFINITIONS "ENABLE_STAT_MONITOR")

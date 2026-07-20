set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(benchmark_numc main.cpp)
target_link_libraries(benchmark_numc PRIVATE numc)

# Enable optimizations for benchmarks
target_compile_options(benchmark_numc PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/O2>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3>
)

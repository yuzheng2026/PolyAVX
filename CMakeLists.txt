cmake_minimum_required(VERSION 3.5)
project(PolyAVX LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_FMA "Enable FMA3 intrinsics" ON)
option(TARGET_ARCH "Target CPU architecture (e.g. native, x86-64)" "native")

if(ENABLE_FMA)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma")
    elseif(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    endif()
endif()

if(TARGET_ARCH)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${TARGET_ARCH}")
    endif()
endif()

add_library(PolyAVX INTERFACE)
target_include_directories(PolyAVX INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    # 示例程序需要链接 cpu_dispatch.cpp
    add_executable(validate validate.cpp cpu_dispatch.cpp)
    target_link_libraries(validate PolyAVX)
    add_executable(bench_basic bench_basic.cpp cpu_dispatch.cpp)
    target_link_libraries(bench_basic PolyAVX)
endif()

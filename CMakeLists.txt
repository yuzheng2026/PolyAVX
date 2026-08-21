cmake_minimum_required(VERSION 3.10)
project(PolyAVX LANGUAGES CXX)

# 添加一个 INTERFACE 库，表示只有头文件，没有可编译源文件
add_library(PolyAVX INTERFACE)

# 让使用 PolyAVX 的目标自动包含当前目录（即 poly_avx.hpp 所在目录）
target_include_directories(PolyAVX INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# 指定 C++98 标准（PolyAVX 兼容 C++98）
target_compile_features(PolyAVX INTERFACE cxx_std_98)

# 可选：为 GCC/Clang 添加优化选项（根据需求调整）
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(PolyAVX INTERFACE -O3 -march=native -mfma)
elseif(MSVC)
    target_compile_options(PolyAVX INTERFACE /O2 /arch:AVX2)
endif()

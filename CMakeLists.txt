
cmake_minimum_required(VERSION 3.16)
project(cfd_lite LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(USE_MPI "Enable MPI support" OFF)
option(USE_CUDA "Enable CUDA support (skeleton only)" OFF)

if (MSVC)
  add_compile_options(/W4 /permissive- /Zc:__cplusplus)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

add_library(core INTERFACE)
target_include_directories(core INTERFACE src)

add_executable(euler2d src/fv_euler2d/main.cpp)
target_link_libraries(euler2d PRIVATE core)

if (USE_MPI)
  find_package(MPI REQUIRED)
  target_compile_definitions(euler2d PRIVATE USE_MPI)
  target_link_libraries(euler2d PRIVATE MPI::MPI_CXX)
endif()

enable_testing()
add_test(NAME run_euler2d_help COMMAND euler2d --help)

cmake_minimum_required (VERSION 2.8.12)
project(extruct CXX)

include (CheckCXXCompilerFlag)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|Intel")
  # Enable C++14
  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)

  if (HAS_CPP14_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  else()
    message(FATAL_ERROR "Unsupported compiler -- C++14 support is needed!")
  endif()
endif()

# Set build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Compile with compiler warnings turned on
if (MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    add_compile_options("/W4")
  endif()
else()
  add_compile_options("-Wall" "-Wextra")
endif()

# Always use libc++ on Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++14")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
endif()

set(PYTHON_EXECUTABLE "/Users/af/anaconda3/bin/python3.6" CACHE FILEPATH "python filepath" FORCE)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/enoki)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tbb)

set(TBB_BUILD_SHARED          ON  CACHE BOOL " " FORCE)
set(TBB_BUILD_STATIC          ON CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC       ON CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC_PROXY ON CACHE BOOL " " FORCE)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/enoki/include
  #${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include
  ${CMAKE_CURRENT_SOURCE_DIR}/tbb/include
  #${CMAKE_CURRENT_SOURCE_DIR}/Eigen
)

#include_directories(pybind11/include)

#add_compile_options("-march=native")
#add_compile_options("-fno-stack-protector")
#add_compile_options("-fomit-frame-pointer")
enoki_set_compile_flags()
enoki_set_native_flags()



add_executable (extruct) 
target_sources(extruct
  PRIVATE 
    main.cpp
    universe.cpp
) 

#pybind11_add_module(extruct main.cpp universe.cpp)

target_link_libraries(extruct PRIVATE tbb fftw3)
#target_link_libraries(1d PRIVATE fftw3)

cmake_minimum_required (VERSION 3.16)

# your project name
project(app)

# include biovoltron CMakeLists.txt
include(${PROJECT_SOURCE_DIR}/Biovoltron/CMakeLists.txt)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include/)

add_definitions( -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" )
# add_compile_options(-isystem ${PROJECT_SOURCE_DIR}/Biovoltron/include/biovoltron/container/xbit_vector.hpp)

# find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
# if(CLANG_TIDY_EXE)
#   set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-p" "${CMAKE_BINARY_DIR}")
#   # set_target_properties(src/main.cpp PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
# endif()

# set(CMAKE_CXX_CLANG_TIDY 'clang-tidy;-extra-arg-before=-config;-extra-arg-before=.clang-tidy')
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_CXX_CLANG_TIDY 
#   clang-tidy;
#   -header-filter=./include/anchor_map.hpp;
#   -checks=*;)


FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
FIND_PACKAGE( Boost 1.73.0 COMPONENTS program_options REQUIRED )
IF (Boost_FOUND)
	# INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}/boost/program_options/)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	message(STATUS ${BOOST_ROOT})
	message(STATUS ${Boost_LIBRARY_DIRS})
	message(STATUS ${BOOST_INCLUDEDIR})
ENDIF()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# link your source code to your project
add_executable(app src/main.cpp)
target_link_libraries(app biovoltron)
# TARGET_LINK_LIBRARIES(app ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(app Boost::program_options)

target_compile_options(app PRIVATE -g)

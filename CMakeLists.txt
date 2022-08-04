set(DEBUG_FLAGS -O0 -g -fno-omit-frame-pointer -pg -fsanitize=address -DDEBUG)
set(PROD_FLAGS -O2 -g -DNDEBUG)
# specify the C++ standard
cmake_minimum_required(VERSION 3.14)

# set the project name
project(MashMap VERSION 3.0)

# set output directory of builds
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set build type
set(CMAKE_BUILD_TYPE RelWithDebInfo)

# set CXX standard
# Things seem to be faster in cxx 20, and there is also std::shift_* 
# Should fall back on 17 if 20 is not supported
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set up linters/checkers
#set(CMAKE_CXX_CPPCHECK cppcheck;--std=c++17;--file-filter=*bitlib*)
#set(CMAKE_CXX_CPPLINT cpplint;--linelength=100;--filter=-whitespace;)
#set(CMAKE_CXX_CLANG_TIDY 
  #clang-tidy;
  #-header-filter=include/;)

# 
# specify global compiler flags
add_compile_options(-march=native -Wfatal-errors -DDEBUG)


# Add targets
file(GLOB MAP_HEADERS "src/map/include" "*.hpp")
file(GLOB MAP_SOURCES "src/map/" "*.cpp")
file(GLOB ALIGN_HEADERS "src/align/include" ".hpp")
file(GLOB ALIGN_SOURCES "src/align/" "*.cpp")
file(GLOB COMMON_HEADERS "src/common/include" ".hpp")
file(GLOB COMMON_SOURCES "src/common/" "*.cpp")
#add_executable(mashmap ${MAP_SOURCES})
add_executable(mashmap-debug "src/map/mash_map.cpp")
add_executable(mashmap "src/map/mash_map.cpp")

################
# Benchmarking #
################
target_compile_options(mashmap PUBLIC ${PROD_FLAGS})

# specify debugging libraries
include_directories(src/)
target_link_libraries(mashmap PRIVATE -lgsl -lgslcblas -lstdc++ -lz -lm -lpthread)

install(TARGETS mashmap DESTINATION .)


###########
# Debugging #
###########
# specify debugging flags
target_compile_options(mashmap-debug PUBLIC ${DEBUG_FLAGS})

# specify debugging libraries
include_directories(src/)
target_link_libraries(mashmap-debug PRIVATE -lgsl -lgslcblas -lstdc++ -lz -lm -lpthread -fsanitize=address)

install(TARGETS mashmap-debug DESTINATION .)




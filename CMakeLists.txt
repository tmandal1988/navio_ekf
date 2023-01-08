# Minimum version of CMake required to build this project
cmake_minimum_required(VERSION 3.16.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME IMUGPSEKF)

# Set name of the project
project(${PROJECT_NAME})


find_package(Threads REQUIRED)

set(PROJECT_SOURCE_DIR
	src
)

set(PROJECT_SOURCE
    ${PROJECT_SOURCE_DIR}/ekf_example_quat.cc
)

set(PROJECT_LIBRARIES_DIR
    libs
)

set(PROJECT_LIBRARIES 
    Matrix
    NavioCommon
    NavioNavio2
    CommonUtils
    Ekf
)

set(PROJECT_LIBRARIES_SUBDIR
    Matrix
    Navio/Common
    Navio/Navio2
    CommonUtils
    Ekf
)

# To enable threading
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)


# Disbale ABI change warning
SET(GCC_NO_ABI_WARNING_FLAG "-Wno-psabi")
# Enable gprof profiling
SET(GCC_TURN_ON_GPROF_FLAG "-pg")
# Set GCC Optimization Level
SET(GCC_OPT_FLAG "-Ofast")
# Set Statis libraries for gprof
# SET(LINKER_ALL_STATIC_FLAG "--static")

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_NO_ABI_WARNING_FLAG} ${GCC_TURN_ON_GPROF_FLAG} ${GCC_OPT_FLAG}")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_NO_ABI_WARNING_FLAG} ${GCC_OPT_FLAG}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LINKER_ALL_STATIC_FLAG}")

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE})

foreach(LIBRARY ${PROJECT_LIBRARIES_SUBDIR})
    add_subdirectory("${PROJECT_LIBRARIES_DIR}/${LIBRARY}")
endforeach(LIBRARY)

target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBRARIES} Threads::Threads)
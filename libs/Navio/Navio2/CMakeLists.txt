cmake_minimum_required(VERSION 3.16.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME
    NavioNavio2Lib
)
set(LIBRARY_NAME
    NavioNavio2
)

file(GLOB LIBRARY_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

file(GLOB LIBRARY_SOURCE CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

project(${PROJECT_NAME})

add_library(${LIBRARY_NAME} STATIC
    ${LIBRARY_HEADERS}
    ${LIBRARY_SOURCE}
)

include_directories(${CMAKE_SOURCE_DIR}/libs/Navio)

target_include_directories(${LIBRARY_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/libs>
  $<INSTALL_INTERFACE:${LIBRARY_NAME}>  # <prefix>/include/${LIBRARY_NAME}
)

target_link_libraries(NavioNavio2 NavioCommon)

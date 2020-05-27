cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME example)

project(${PROJECT_NAME})

if (WIN32)

    set(LIB_NAME libVectorization_parallelograms.dll)

elseif (UNIX)

    set(LIB_NAME libVectorization_parallelograms.so)

else ()

    message(FATAL_ERROR "Unknown operating system")

endif ()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../include/)
set(LIBS_DIR ${PROJECT_SOURCE_DIR}/../)

include_directories(${INCLUDE_DIR})
link_directories(${LIBS_DIR})

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${LIB_NAME})

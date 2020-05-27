cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME Vecorization_parallelograms)

project(${PROJECT_NAME})

set(SOURCES_DIR ${PROJECT_SOURCE_DIR}/src/main/)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/)

if (WIN32)

    set(OPENCV_LIBS_DIR ${PROJECT_SOURCE_DIR}/lib/windows)
    set(OPENCV_LIBS opencv_world420.lib opencv_world420d.lib)

elseif (UNIX)

    set(OPENCV_LIBS_DIR ${PROJECT_SOURCE_DIR}/lib/linux)
    set(OPENCV_LIBS libopencv_core.so libopencv_imgproc.so)

else ()

    message(FATAL_ERROR "Unknown operating system")

endif ()

include_directories(${INCLUDE_DIR})
link_directories(${OPENCV_LIBS_DIR})

set(SOURCES ${SOURCES_DIR}main.cpp
        ${SOURCES_DIR}LinesDetection.cpp
        ${SOURCES_DIR}Vectorization_parallelograms.cpp
        ${SOURCES_DIR}MSD.cpp)

set(HEADERS ${INCLUDE_DIR}debug.h
        ${INCLUDE_DIR}Vectorization_parallelograms.h
        ${INCLUDE_DIR}MSD.h
        ${INCLUDE_DIR}LinesDetection.h
        ${INCLUDE_DIR}structures.h)

add_executable(${PROJECT_NAME} ${HEADER} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${OPENCV_LIBS})

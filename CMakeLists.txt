cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME Vecorization_parallelograms)

project(${PROJECT_NAME})

if (WIN32)
    set(SOURCES_DIR ${PROJECT_SOURCE_DIR}\\src\\main\\)
elseif (UNIX)
    set(SOURCES_DIR ${PROJECT_SOURCE_DIR}/src/main/)
else ()
    message(FATAL_ERROR "����������� �������")
endif ()

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SOURCES_DIR}LinesDetection.cpp
        ${SOURCES_DIR}Vectorization_parallelograms.cpp
        ${SOURCES_DIR}MSD.cpp
        ${SOURCES_DIR}main.cpp
        ${SOURCES_DIR}debug.h
        ${SOURCES_DIR}Vectorization_parallelograms.h
        ${SOURCES_DIR}MSD.h
        ${SOURCES_DIR}LinesDetection.h
        ${SOURCES_DIR}Structures.h)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
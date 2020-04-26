cmake_minimum_required(VERSION 2.8)
set(PROJECT_NAME Vecorization_parallelograms)

project(${PROJECT_NAME})

if (WIN32)
    set(SOURCES_DIR ${PROJECT_SOURCE_DIR}\\src\\main\\)
elseif (UNIX)
    set(SOURCES_DIR ${PROJECT_SOURCE_DIR}/src/main/)
else ()
    message(FATAL_ERROR "Неизвестная система")
endif ()

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SOURCES_DIR}Vectorization_OpenCV.cpp
        ${SOURCES_DIR}Vectorization_parallelograms.cpp
        ${SOURCES_DIR}MSD.cpp
        ${SOURCES_DIR}main.cpp src/main/debug.h)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

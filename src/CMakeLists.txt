cmake_minimum_required(VERSION 4.1)
project(NS)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
        cli11
        URL https://github.com/CLIUtils/CLI11/releases/download/v2.4.1/CLI11.hpp
        DOWNLOAD_NO_EXTRACT TRUE
)

FetchContent_MakeAvailable(cli11)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
)
list(FILTER SOURCES EXCLUDE REGEX "/cmake-build-[^/]*/")


set(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
foreach(src IN LISTS SOURCES)
    get_filename_component(dir "${src}" DIRECTORY)
    list(APPEND INCLUDE_DIRS "${dir}")

endforeach()
list(REMOVE_DUPLICATES INCLUDE_DIRS)

add_executable(NS ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
        ${INCLUDE_DIRS}
        ${cli11_SOURCE_DIR}
)
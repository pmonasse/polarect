cmake_minimum_required(VERSION 2.8)

project(Polarect C CXX)

enable_testing()
add_subdirectory(third_party/CppUnitLite)
include_directories(third_party)
macro (UNIT_TEST NAME EXTRA_LIBS)
  add_executable(${NAME}_test ${NAME}_test.cpp)
  target_link_libraries(${NAME}_test
                        ${EXTRA_LIBS} # Extra libs MUST be first.
                        CppUnitLite)
  add_test(NAME ${NAME} COMMAND ${NAME}_test)
endmacro()

add_subdirectory(libOrsa)
add_subdirectory(libSplinter)
add_subdirectory(libImageIO)
add_subdirectory(third_party/sift_anatomy)

# Executable demo to estimate the fundamental matrix from SIFT matches
add_executable(polarect
               polarect.cpp polarect.h
               main.cpp
               siftMatch.hpp)
target_link_libraries(polarect Orsa Splinter ImageIO SiftAnatomy)

if(UNIX)
    set_target_properties(polarect PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra")
endif()

cmake_minimum_required(VERSION 3.20)

project(solenoid
   VERSION 0.1.0
   DESCRIPTION "Tabajo en grupo: solenoide y permeabilidad magnética"
   LANGUAGES CXX
)

add_executable(${PROJECT_NAME} main.cpp)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_23)

set_target_properties(${PROJECT_NAME} PROPERTIES
   CXX_EXTENSIONS OFF
   RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug
   RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release
)

#-------------------------------------------------------------------
# integración de bibliotecas de terceros:

set (gsl-lite_DIR "C:/Program Files (x86)/gsl_lite/lib/cmake/gsl-lite")

find_package(Boost REQUIRED)
find_package(cppexten REQUIRED)
find_package(dlib REQUIRED)
find_package(Matplot++ REQUIRED)
find_package(mp-units REQUIRED)

target_link_libraries(${PROJECT_NAME} PRIVATE
   Boost::headers
   cppexten::cppexten
   dlib::dlib
   Matplot++::matplot
   mp-units::mp-units
)

#-------------------------------------------------------------------
# política de avisos específica de GCC:

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   target_compile_options(${PROJECT_NAME} PRIVATE /W3 /WX)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
   target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Werror)
endif()

#-------------------------------------------------------------------
# mensajes a emitir en construcción:

message("-- CMake Generator: ${CMAKE_GENERATOR}")

get_target_property(CFEATURES ${PROJECT_NAME} COMPILE_FEATURES)
message("-- Target compile features: ${CFEATURES}")

get_target_property(COPTIONS ${PROJECT_NAME} COMPILE_OPTIONS)
message("-- Target compile options: ${COPTIONS}")
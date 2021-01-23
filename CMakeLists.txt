cmake_minimum_required(VERSION 3.11)
project(Infrared)

set(SOURCE_DIR "src/infrared")
include_directories(
      ${SOURCE_DIR}
)



# include(pybind11.cmake)
find_package(pybind11 REQUIRED)
pybind11_add_module(libinfrared ${SOURCES} "${SOURCE_DIR}/libinfrared.cpp")

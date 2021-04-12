cmake_minimum_required(VERSION 3.14)

# project(spdlog_cmake)
project(MyTests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
# for experimental file system
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -std=c++17 -lstdc++fs")

add_subdirectory(spdlog)
# add_subdirectory(myLearn)

# find_package(spdlog)

add_executable(simpleLog myLearn/simpleLogger.cpp)

include_directories(include)

set(HEALPIX_SOURCE_DIR "/home/panwei/Documents/Healpix_3.50/src/cxx/basic_gcc/")
set(CFITSIO_SOURCE_DIR "/home/panwei/Documents/cfitsio-3.47/")
include_directories(${HEALPIX_SOURCE_DIR}/inc)
link_directories(${HEALPIX_SOURCE_DIR}/lib ${CFITSIO_SOURCE_DIR}/lib)


add_subdirectory(source)

set(CMAKE_CURRENT_SOURCE_DIR "./source")
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} MAIN_SOURCE)

# add_library(points ./sphericalPoint.cpp ./)
add_executable(test test.cpp)
target_link_libraries(test sphere2pCorr fileIO)
target_link_libraries(test stdc++fs)
target_link_libraries(test spdlog)


add_executable(exampleTest exampleTest.cpp)
target_link_libraries(exampleTest sphere2pCorr fileIO testFunc)
target_link_libraries(exampleTest stdc++fs)
target_link_libraries(exampleTest spdlog)

add_executable(exampleIOtest exampleIOtest.cpp)
target_link_libraries(exampleIOtest testFunc)

target_link_libraries(simpleLog
  spdlog
)

link_directories(${CMAKE_CURRENT_LIST_DIR}/lib)

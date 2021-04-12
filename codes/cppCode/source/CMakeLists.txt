include_directories(${CMAKE_SOURCE_DIR}/include)  # include all the headers
include_directories(${CMAKE_SOURCE_DIR}/spdlog/include)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")


set(HEALPIX_SOURCE_DIR "/home/panwei/Documents/Healpix_3.50/src/cxx/basic_gcc/")
set(CFITSIO_SOURCE_DIR "/home/panwei/Documents/cfitsio-3.47/")
include_directories(${HEALPIX_SOURCE_DIR}/inc)
link_directories(${HEALPIX_SOURCE_DIR}/lib ${CFITSIO_SOURCE_DIR}/lib)
# link_directories(${CFITSIO_SOURCE_DIR}/lib)
set(LIBRARY_SOURCE "sharp  healpix_cxx cxxsupport c_utils  pocketfft cfitsio") #  healpix_cxx cxxsupport c_utils  pocketfft cfitsio
# file IO library
add_library(fileIO fileIO.cpp)
add_library(points poissonPoint2d.cpp sphericalPoint.cpp)
add_library(HealpixCorr healpixCorrelation.cpp binsGenerator.cpp)
add_library(sphere2pCorr sphereTwoCorrelation.cpp findPair.cpp)
add_library(testFunc testFunc.cpp)

# add library dependencies

# add pthread dependencies
find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(points PUBLIC "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
  target_link_libraries(points "${CMAKE_THREAD_LIBS_INIT}")
endif()
message(${LIBRARY_SOURCE})
target_link_libraries(points spdlog)
target_link_libraries(HealpixCorr points sharp healpix_cxx cxxsupport c_utils pocketfft cfitsio)
# target_link_libraries(HealpixCorr points ${LIBRARY_SOURCE}) # need to include the source library
target_link_libraries(sphere2pCorr HealpixCorr)
# target_link_libraries(fileIO HealpixCorr)
target_link_libraries(testFunc HealpixCorr fileIO)




cmake_minimum_required(VERSION 3.3)

project(HDFacet)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)
find_package(Matlab)
find_package(Threads)

if(APPLE)
	find_library(LIB_MATLAB_DATA_ARRAY 
		NAMES MatlabDataArray
		PATHS ${Matlab_ROOT_DIR}/extern/bin/maci64
		NO_DEFAULT_PATH)

	find_library(LIB_MATLAB_ENGINE
		NAMES MatlabEngine
		PATHS ${Matlab_ROOT_DIR}/extern/bin/maci64
        NO_DEFAULT_PATH)
else()
	message(FATAL_ERROR "System not MacOS, help cmake find Matlab engine libraries for C++.")
endif()

include(${CGAL_USE_FILE})
include_directories(${Matlab_INCLUDE_DIRS})


include_directories(.)
include_directories(eigen)

file(GLOB Headers "src/*.h")
file(GLOB Sources "src/*.cpp")

add_executable(HDFacet ${Headers} ${Sources})

target_link_libraries( HDFacet 
	${CGAL_LIBRARIES} 
	${CGAL_3RD_PARTY_LIBRARIES} 
	${LIB_MATLAB_DATA_ARRAY}
	${LIB_MATLAB_ENGINE}
	Threads::Threads 	
)	

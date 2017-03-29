# -*- CMakeLists.txt. Do not edit by hand -*-
# 　Copyright 2013 by Hao Wang  
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

cmake_minimum_required(VERSION 2.8.11)

project(Hieralign)

# Define some variables
set(PROJECT_PATH "${CMAKE_SOURCE_DIR}/Hieralign")


include_directories(
    .
    .

)


# Compiler options
add_definitions(-W)
add_definitions(-fopenmp)
add_definitions(-std=c++11) 
add_definitions(-Wall)
add_definitions(-O3)
add_definitions(
    -DNDEBUG
)


# Linker options
set(LINK_OPTIONS -fopenmp)


if(WIN32)
    # Resource options
endif(WIN32)

# Library path
set(CMAKE_LDFLAGS "${CMAKE_LDFLAGS} -L. ")

# Define the CXX sources
set ( CXX_SRCS
    ${CMAKE_SOURCE_DIR}/src/main.cc
    ${CMAKE_SOURCE_DIR}/src/IBM1.cc
    ${CMAKE_SOURCE_DIR}/src/Hieralign.cc
    ${CMAKE_SOURCE_DIR}/src/VBIBM1.cc
)

set_source_files_properties(
    ${CXX_SRCS} PROPERTIES COMPILE_FLAGS 
    " -W -fopenmp -std=c++11 -O3 -fopenmp -std=c++11")

if(WIN32)
    enable_language(RC)
    set(CMAKE_RC_COMPILE_OBJECT
        "<CMAKE_RC_COMPILER> ${RC_OPTIONS} -O coff -i <SOURCE> -o <OBJECT>")
endif(WIN32)

add_executable(Hieralign ${RC_SRCS} ${CXX_SRCS} ${C_SRCS})
target_link_libraries(Hieralign ${LINK_OPTIONS})



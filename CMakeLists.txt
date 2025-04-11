cmake_minimum_required(VERSION 3.10)
project(HOOMDClusterParser CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBXML2 REQUIRED libxml-2.0)

include_directories(${LIBXML2_INCLUDE_DIRS})
include_directories(include)

add_executable(hoomd_cluster2
    src/main.cpp
    src/parser.cpp
    src/clustering.cpp
)

add_executable(clout_ana2
    tools/clout_ana.cpp
    tools/analyse_clfiles.cpp
)

target_link_libraries(hoomd_cluster2 ${LIBXML2_LIBRARIES} m)

install(TARGETS hoomd_cluster2 DESTINATION bin)
install(TARGETS clout_ana2 DESTINATION bin)
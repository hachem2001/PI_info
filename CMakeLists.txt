cmake_minimum_required(VERSION 3.17)

project(PI_info VERSION 1.0 DESCRIPTION "PI Info." LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set_property(GLOBAL PROPERTY USE_FOLDERS ON) # Useful for IDEs.

# Source code here. Includes both reusable source code, and the application of it.
add_subdirectory(src)



#
# Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
#
# $Id: CMakeLists.txt 163 2012-06-26 14:22:38Z hhl $
#



cmake_minimum_required (VERSION 2.8)
set (CMAKE_LEGACY_CYGWIN_WIN32 0)

if (NOT UNIX)
  message (FATAL_ERROR "Non Unix-like environments currently not supported")
endif (NOT UNIX)

project (molprep C Fortran)

set (VERSION_MAJOR 0)
set (VERSION_MINOR 2)
set (VERSION_PATCH 2)

add_subdirectory(src)
add_subdirectory(data)

set (CMAKE_INSTALL_PREFIX "/usr/local")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Options are: None Debug Release."
       FORCE)
endif (NOT CMAKE_BUILD_TYPE)

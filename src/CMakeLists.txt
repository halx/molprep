#
# Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
#
# $Id: CMakeLists.txt 166 2012-07-05 08:46:10Z hhl $
#



set (EXTRA_LIBS "")

find_package(ZLIB)

if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  set (EXTRA_LIBS ${EXTRA_LIBS} ${ZLIB_LIBRARIES})
  set (HAVE_ZLIB 1)
endif (ZLIB_FOUND)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
)

include_directories(${PROJECT_BINARY_DIR})

add_executable(molprep molprep.c hbuild.c pdb.c protonate.c ssbuild.c top.c
               propka/propka.F)

target_link_libraries(molprep ${EXTRA_LIBS} molprep_util)

install (TARGETS molprep DESTINATION bin)

# Open64 is detected as GNU too
if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
  # pcc is recognized as GNU compiler
  if (CMAKE_C_COMPILER MATCHES "pcc")
    set (CMAKE_C_FLAGS "")
    set (CMAKE_C_FLAGS_DEBUG "-g3")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
  else ()
    # we assume that Clang accepts the same flags as GCC
    set (CMAKE_C_FLAGS "-std=c99 -pedantic")
    set (CMAKE_C_FLAGS_DEBUG "-g3 -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fno-common -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fno-common -std=c99 -pedantic")
    set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
  endif (CMAKE_C_COMPILER MATCHES "pcc")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_C_FLAGS "-std=c99 -strict-ansi")
  set (CMAKE_C_FLAGS_DEBUG "-g3 -w2 -Wbrief -Wall -Wcheck -check-uninit")
  set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")
  set (CMAKE_C_FLAGS "-Xa -c99")
  set (CMAKE_C_FLAGS_DEBUG "-g -Mbounds -Minform=warn")
# optimisation crashes the C code on the 64bit platforms
  set (CMAKE_C_FLAGS_RELEASE "-g -DNDEBUG")
elseif (CMAKE_C_COMPILER MATCHES "pcc")
  # pcc is recognized as GNU so this will never be checked
  # but pcc appears to have problem with macro expansion
  set (CMAKE_C_FLAGS "")
  set (CMAKE_C_FLAGS_DEBUG "-g3")
  set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
elseif (CMAKE_C_COMPILER_ID STREQUAL "PathScale")
  set (CMAKE_C_FLAGS "-std=c99 -pedantic")
  set (CMAKE_C_FLAGS_DEBUG "-g3")
  set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Cray")
  set (CMAKE_C_FLAGS "-h c99")
  set (CMAKE_C_FLAGS_DEBUG "-Gf")
  set (CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG -h msglevel_4 -O2")
elseif (CMAKE_C_COMPILER_ID STREQUAL "SunPro")
  set (CMAKE_C_FLAGS "-xc99")
  set (CMAKE_C_FLAGS_DEBUG "-g")
  set (CMAKE_C_FLAGS_RELEASE "-xO2 -DNDEBUG")
else ()
  message (WARNING "Unsupported compiler: ${CMAKE_C_COMPILER}")
endif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID STREQUAL "Clang")


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Open64 4.2.4 doesn't know about '-finit-local-zero'
  set (CMAKE_Fortran_FLAGS "-finit-local-zero")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g3 -Wall -Wextra -Wcharacter-truncation -Wtabs -Wunderflow -Wshadow -Wcast-align -Waggregate-return -fno-common -Wcharacter-truncation -Wconversion-extra -Wunderflow -Wunused-parameter -fcheck=all")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -Wno-unused -Wno-unused-parameter")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
# either add "-heap-arrays [n]" or increase stack size via "ulimit"
  set (CMAKE_Fortran_FLAGS "")
# "-check all" crashes with a segfault
  set (CMAKE_Fortran_FLAGS_DEBUG "-g3")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_LINK_FLAGS "-nofor-main")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set (CMAKE_Fortran_FLAGS "")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -Mbounds -Minform=warn")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_LINK_FLAGS "-Mnomain")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PathScale")
# either add '-LANG:heap_allocation_threshold=10240' or increase stack size via "ulimit"
  set (CMAKE_Fortran_FLAGS "-zerouv")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g3")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
# set_target_properties(molprep PROPERTIES LINKER_LANGUAGE C)
  set (CMAKE_Fortran_FLAGS "-e 0 -h zero")
  set (CMAKE_Fortran_FLAGS_DEBUG "-G2")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -m 4")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
# either add 'stackvar=no' which is the default or increase stack size via "ulimit"
  set (CMAKE_Fortran_FLAGS "-xcheck=init_local")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
else ()
  message (WARNING "Unsupported compiler: ${CMAKE_Fortran_COMPILER}")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(util)

# interface checking fails for Cray and Open64
# with Cray compiler comment the following lines out, use
# 'set_target_properties' above and remove '-lomp_lib' from link.txt
include (FortranCInterface)
FortranCInterface_VERIFY()

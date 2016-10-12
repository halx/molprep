#
# Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
#
# $Id: CMakeLists.txt 162 2012-06-25 14:33:29Z hhl $
#



add_library(molprep_util STATIC llist.c darray.c hashtab.c hashfuncs.c util.c
                                zio.c)

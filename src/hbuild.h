/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: hbuild.h 161 2012-06-25 12:51:40Z hhl $
 *
 */


#ifndef _HBUILD_H
#define _HBUILD_H      1

#include "util/hashtab.h"

void hbuild(pdb_root *pdb, const Hashtable *top, char altLoc);

#endif

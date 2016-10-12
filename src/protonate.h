/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: protonate.h 159 2012-06-25 08:42:40Z hhl $
 *
 */


#ifndef _PROTONATE_H
#define _PROTONATE_H      1

void protonate(pdb_root *pdb, const Hashtable *top, const char *ttb_filename,
	       float pH, char altLoc);

#endif

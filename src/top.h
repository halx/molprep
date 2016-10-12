/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: top.h 161 2012-06-25 12:51:40Z hhl $
 *
 */


#ifndef _TOP_H
#define _TOP_H      1

#include "pdb.h"
#include "util/hashtab.h"


typedef struct _topol_hydro {
  unsigned int nhyd;
  unsigned int type;
  float xhdist;
  char atoms[5][PDB_ATOM_NAME_LEN];
} topol_hydro;

typedef struct _topol {
  char res_type;
  char mol_type;
  char resName[PDB_RES_NAME_LEN];
  struct _topol *first_term;
  struct _topol *last_term;
  char **heavy_atoms;
  topol_hydro **hydrogens;
} topol;

typedef struct _topol_hash {
  topol *data;
  Hashtable *hash_table;
} topol_hash;

int topcmp(const void *p1, const void *p2);
topol_hash *top_read(topol_hash* top_hash, const char *filename);
void top_destroy(topol_hash *top);

#ifndef NDEBUG
void top_print(topol *top);
#endif

#endif

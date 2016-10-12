/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: pdb.h 161 2012-06-25 12:51:40Z hhl $
 *
 */



#ifndef _PDB_H
#define _PDB_H      1

#include "util/vec.h"

#define PDB_LINE_LEN 82
#define PDB_RES_NAME_LEN 5
#define PDB_ATOM_NAME_LEN 5
#define PDB_SEG_NAME_LEN 5
#define PDB_SERIAL_LEN 6
#define PDB_ELEMENT_LEN 3
#define PDB_CHARGE_LEN 3
#define PDB_ID_LEN 5
#define PDB_SSBOND_SYMOP_LEN 7


struct _ssbond {
  char chainID;
  char icode;
  int seqNum;
  char SymOP[PDB_SSBOND_SYMOP_LEN];
};

typedef struct _pdb_ssbond {
  int serNum;
  float Length;
  struct _ssbond ss1;
  struct _ssbond ss2;
} pdb_ssbond;

typedef struct _pdb_atom {
  char serial[PDB_SERIAL_LEN];	/* actually int but unreliable */
  char name[PDB_ATOM_NAME_LEN];
  char altLoc;
  fvec pos;
  float occupancy;
  float tempFactor;
  char element[3];
  char charge[3];
  struct _pdb_residue *residue;
  struct _pdb_atom *next;
} pdb_atom;

typedef struct _pdb_residue {
  char iCode;
  int resSeq;
  char rectype;
  char resName[PDB_RES_NAME_LEN];
  char segID[PDB_SEG_NAME_LEN];	/* old PDB v2.2, CHARMM still uses it */
  pdb_atom *first_atom;
  struct _pdb_chain *chain;
  struct _pdb_residue *next;
} pdb_residue;

typedef struct _pdb_chain {
  char chainID;
  pdb_residue *first_residue;
  struct _pdb_chain *next;
} pdb_chain;

typedef struct _pdb_root {
  unsigned int natoms;
  unsigned int nres;
  unsigned int nchains;
  int model_no;
  char ID[PDB_ID_LEN];
  char cryst1[PDB_LINE_LEN];
  pdb_ssbond **ssbonds;
  pdb_chain *first_chain;
} pdb_root;


pdb_atom *pdb_insert_atom_node(pdb_atom *old);

int pdb_format_atom(char *restrict dest, const char *restrict src);
char *pdb_format_residue(char *restrict dest, const char *restrict src);

pdb_root *pdb_read(pdb_root *pdb, const char *filename,  const char* ss_name,
		   int model_no, int *nssb);
void pdb_write(pdb_root *pdb, const char *filename, const char *format,
	       const char* ss_name, char altLoc);
void pdb_destroy(pdb_root *pdb);

#endif

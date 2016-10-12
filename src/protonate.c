/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Protonate protein via interface to PROPKA 2.0 (in Fortran 77)
 *
 * VERY IMPORTANT: PROPKA is ABSOLUTELY RELIANT on a 'well-behaved' PDB file:
 * complete side chains; atoms in 'PDB order'; only "ATOM", "HETATM", and "LG"
 * records
 *
 *
 * $Id: protonate.c 163 2012-06-26 14:22:38Z hhl $
 *
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "common.h"
#include "pdb.h"
#include "top.h"
#include "util/hashtab.h"
#include "util/util.h"
#include "propka/propka.h"


#define PROPKA_PDB_FILE "propka.pdb"
#define PROPKA_OUT_FILE "propka.out"
#define RETURN_STRING_SIZE 1048576 // possibly good for up to ~800.000 heavy atoms
#define PDB_PROPKA_FORMAT "%-6s%5s %4s%c%4s%c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n"

#define TTB_DELIMITER  " =->\t\n"
#define TTB_LINE_LEN 82
#define AR_TAB_SIZE 12	      // atoms per residue table
#define TITR_TAB_SIZE 7	      // the number of titratable residues as below
#define TITRATABLE(res) (STREQ(res, "ARG ") || STREQ(res, "ASP ") || \
			 STREQ(res, "CYS ") || STREQ(res, "GLU ") || \
			 STREQ(res, "HIS ") || STREQ(res, "LYS ") || \
			 STREQ(res, "TYR ") )

struct _titr_table {
  char name[PDB_ATOM_NAME_LEN];
  char prot_name[PDB_ATOM_NAME_LEN];
};

struct _propka_table {
  char chainID;
  int resSeq;
  float pKa;
  char resName[PDB_RES_NAME_LEN];
  char *prot_name;
};


void protonate(pdb_root *pdb, const Hashtable *top, const char *ttb_filename,
	       float pH, char altLoc)
{
  int resSeq = 0, retc;
  unsigned int i, maxar;
  unsigned int serno = 0;
  unsigned int atom_cnt = 0, residue_cnt = 0, titr_cnt = 0, line_cnt = 0;
  unsigned int nfields, nttb;

  bool found;

  char chainID = ' ';
  char resName[PDB_RES_NAME_LEN];
  char serial[6];
  char pdb_file[] = PROPKA_PDB_FILE;
  char out_file[] = PROPKA_OUT_FILE;
  char tmp[PDB_ATOM_NAME_LEN];
  char buffer[TTB_LINE_LEN];
  char return_string[RETURN_STRING_SIZE];  // NOTE: do NOT allocate from heap!
  char *key, *val, *bufp, **heavy;

  float pKa = 0.0;

  pdb_atom *curr_atom;
  pdb_residue *curr_residue;
  pdb_chain *curr_chain;

  topol *top_entry;

  Hashnode *curr_node;

  struct _propka_table *propka_table, *tp;

  // each atom of the following residues is stored individually:
  // ASP, GLU, ARG, CYS, HIS, LYS, TYR, GLN, ASN, TRP, SER, THR
  unsigned int ar[AR_TAB_SIZE] = {0};

  struct _titr_table ttb[TITR_TAB_SIZE] = {
    {"\0", "\0"}, {"\0", "\0"}, {"\0", "\0"}, {"\0", "\0"}, {"\0", "\0"},
    {"\0", "\0"}, {"\0", "\0"}};

  FILE *pdb_stream, *ttb_stream;



  if (!(pdb_stream = fopen(PROPKA_PDB_FILE, "w")) ) {
    perror(PROPKA_PDB_FILE);
    exit(2);
  }


  for (curr_chain = pdb->first_chain; curr_chain; curr_chain = curr_chain->next) {
    for (curr_residue = curr_chain->first_residue;
	 curr_residue && curr_residue->chain == curr_chain;
	 curr_residue = curr_residue->next) {  // residue


      if ( !(curr_node = hash_search(top, curr_residue->resName,
				     strlen(curr_residue->resName) ) ) ) {
	continue;
      }

      top_entry = hash_node_get_data(curr_node);

      // skip non-protein residues
      if (top_entry->mol_type != 'P' ||
	  curr_residue->rectype != 'A') {
	continue;
      }

      // FIXME: C-terminus may need OXT
      for (heavy = top_entry->heavy_atoms; *heavy; heavy++) {  // heavy
	found = false;

	for (curr_atom = curr_residue->first_atom;
	     curr_atom && curr_atom->residue == curr_residue;
	     curr_atom = curr_atom->next) {

	  if (curr_atom->altLoc != altLoc && curr_atom->altLoc != ' ') {
	    continue;
	  }

	  if (STRNEQ(*heavy, curr_atom->name, PDB_ATOM_NAME_LEN-1) ) {
	    found = true;
	    break;
	  }
	}

	if (found) {
	  atom_cnt++;
	  serno++;

	  if (serno > 99999)
	    serno = 1;

	  sprintf(serial, "%i", serno);

	  fprintf(pdb_stream, PDB_PROPKA_FORMAT, "ATOM",
		  serial, curr_atom->name, curr_atom->altLoc,
		  curr_residue->resName, curr_chain->chainID,
		  curr_residue->resSeq, curr_residue->iCode,
		  curr_atom->pos[0], curr_atom->pos[1], curr_atom->pos[2],
		  curr_atom->occupancy, curr_atom->tempFactor);
	} else {
	  prwarn("PROPKA cannot protonate: incomplete amino acid (%s %d%c %c)\n",
		 curr_residue->resName, curr_residue->resSeq, curr_residue->iCode,
		 curr_chain->chainID);

	  fclose(pdb_stream);
	  unlink(PROPKA_PDB_FILE);

	  return;
	}
      }	// heavy

      // number of atoms per residue stored individually
      if (STREQ(curr_residue->resName, "ASP ") ) {
	ar[0]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "GLU ") ) {
	ar[1]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "ARG ") ) {
	ar[2]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "CYS ") ) {
	ar[3]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "HIS ") ) {
	ar[4]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "LYS ") ) {
	ar[5]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "TYR ") ) {
	ar[6]++;
	titr_cnt++;
      } else if (STREQ(curr_residue->resName, "GLN ") ) {
	ar[7]++;
      } else if (STREQ(curr_residue->resName, "ASN ") ) {
	ar[8]++;
      } else if (STREQ(curr_residue->resName, "TRP ") ) {
	ar[9]++;
      } else if (STREQ(curr_residue->resName, "SER ") ) {
	ar[10]++;
      } else if (STREQ(curr_residue->resName, "THR ") ) {
	ar[11]++;
      }

      residue_cnt++;
    } // residue
  } // chain

  fclose(pdb_stream);

  if (atom_cnt < 5) {		// GLY is smallest amino acid with 4 bb atoms
    prwarn("PROPKA cannot protonate: PDB does not contain protein/peptide\n");

    unlink(PROPKA_PDB_FILE);

    return;
  }

  unlink(PROPKA_OUT_FILE);	// PROPKA requires 'NEW' status

  if (titr_cnt < 1) {
    unlink(PROPKA_PDB_FILE);
    return;
  }

  prnote("PROPKA 2.0 will analyze %i titratable residues (out of %i) at pH %.2f\n",
	 titr_cnt, residue_cnt, pH);

  titr_cnt += 3;		// extra space to accomodate "N+", "C-"

  if (16 * titr_cnt >= RETURN_STRING_SIZE) { // 16 characters/site from propka.F
    prerror(4, "Current protonation table too small (%i). "
	    "Number of titratable sites is %i\n", RETURN_STRING_SIZE,
	    titr_cnt - 3);
  }

  maxar = ar[0] + ar[1];	// PROPKA stores ASP and GLU together

  for (unsigned int idx = 2; idx < AR_TAB_SIZE; idx++) {
    if (ar[idx] > maxar) {
      maxar = ar[idx];
    }
  }

  // ASP+GLU may store OXT atom, LYS may store terminal N atom
  maxar++;


  // At this point we must have a 'clean' PDB ready for PROPKA

  retc = runpka_(&atom_cnt, &residue_cnt, &maxar,
		 pdb_file, out_file, return_string,
		 strlen(pdb_file), strlen(out_file), RETURN_STRING_SIZE);

  if (retc) {
    prwarn("PROPKA cannot protonate.\n");
    unlink(PROPKA_PDB_FILE);
    return;
  }


  prnote("reading titratable translation table from %s\n", ttb_filename);

  if (!(ttb_stream = fopen(ttb_filename, "r")) ) {
    perror(ttb_filename);
    exit(EXIT_FAILURE);
  }

  i = 0;

  while (fgets(buffer, TTB_LINE_LEN, ttb_stream) ) {  // read lines
    line_cnt++;

    if (!(bufp = normln(buffer)) )
      continue;

    if ( !(key = strtok(bufp, TTB_DELIMITER) ) ) {
      prerror(4, "%s: no key found in input (line %d).\n",
	      ttb_filename, line_cnt);
    }

    if ( !(val = strtok(NULL, TTB_DELIMITER) ) ) {
      prerror(4, "%s: no value found in input (line %d).\n", 
	      ttb_filename, line_cnt);
    }

    if (!pdb_format_residue(tmp, key) ) {
      prerror (2, "%s: atom name %s too long in line %d.\n", ttb_filename, key,
	       line_cnt);
    }

    if (!TITRATABLE(tmp) ) {
      prwarn ("%s: %s is not a titrable site in line %d.\n", ttb_filename, tmp,
	      line_cnt);
      continue;
    }

    strncpy(ttb[i].name, tmp, PDB_ATOM_NAME_LEN-1);
    ttb[i].name[PDB_ATOM_NAME_LEN-1] = '\0';

    if (!pdb_format_residue(ttb[i].prot_name, val) ) {
      prerror (2, "%s: atom name %s too long in line %d.\n", ttb_filename, val,
	       line_cnt);
    }

    i++;
  }

  fclose(ttb_stream);

  nttb = i;

  propka_table = allocate(titr_cnt * sizeof *propka_table);
  bufp = return_string;

  i = 0;

  while ( (bufp = strtok(bufp, "|") ) ) {
    nfields = sscanf(bufp, "%3s%4i%c%f", resName, &resSeq, &chainID, &pKa);

    if (nfields < 4) {
      prerror(3, "error in parsing propka output\n");
    }

    if (!pdb_format_residue(propka_table[i].resName, resName) ) {
      prerror(3, "propka output: atom name %s too long.\n", resName);
    }

    assert(i < titr_cnt);

    propka_table[i].resSeq = resSeq;
    propka_table[i].chainID = chainID;
    propka_table[i].pKa = pKa;
    propka_table[i].prot_name = NULL;

    for (unsigned int j = 0; j < nttb; j++) {
      if (STRNEQ(propka_table[i].resName, ttb[j].name, PDB_ATOM_NAME_LEN-1)) {
	propka_table[i].prot_name = ttb[j].prot_name;
      }
    }

    i++;
    bufp = NULL;
  }

  propka_table[i].resName[0] = '\0';


  for (curr_chain = pdb->first_chain; curr_chain; curr_chain = curr_chain->next) {
    for (curr_residue = curr_chain->first_residue;
	 curr_residue && curr_residue->chain == curr_chain;
	 curr_residue = curr_residue->next) {  // residue

      if (!TITRATABLE(curr_residue->resName) )
	continue;

      for (tp = propka_table; *(tp->resName); tp++) {
	if (tp->prot_name && STREQ(tp->resName, curr_residue->resName) &&
	    tp->resSeq == curr_residue->resSeq &&
	    tp->chainID == curr_chain->chainID) {

	  // FIXME: how to deal with termini? "N+" and "C-" ignored at the moment
	  if (pH < tp->pKa) {
	    if (STREQ(curr_residue->resName, "HIS ") ||
		STREQ(curr_residue->resName, "ASP ") ||
		STREQ(curr_residue->resName, "GLU ") ) {

	      prnote("protonating %s %i %c (pKa = %.2f)\n",
		     curr_residue->resName, curr_residue->resSeq,
		     curr_chain->chainID, tp->pKa);

	      strncpy(curr_residue->resName, tp->prot_name, PDB_RES_NAME_LEN-1);
	      curr_residue->resName[PDB_RES_NAME_LEN-1] = '\0';
	    }
	  } else {
	    if (STREQ(curr_residue->resName, "LYS ") ||
		STREQ(curr_residue->resName, "CYS ") || 
		STREQ(curr_residue->resName, "ARG ") ||
		STREQ(curr_residue->resName, "TYR ") ) {

	      prnote("deprotonating %s %i %c (pKa = %.2f)\n",
		     curr_residue->resName, curr_residue->resSeq,
		     curr_chain->chainID, tp->pKa);

	      strncpy(curr_residue->resName, tp->prot_name, PDB_RES_NAME_LEN-1);
	      curr_residue->resName[PDB_RES_NAME_LEN-1] = '\0';
	    }
	  }

	  break;
	}
      }
    } // residue
  } // chain

  free(propka_table);
}

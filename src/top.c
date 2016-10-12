/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Read a topology database file.  Molecule types are tagged as [str] where str
 * is one of 'proteins', 'DNA', 'RNA' or 'other'.  Residue names followed by
 * optional alias names are preceded with RESIDUE.  FTERM and LTERM link to the
 * first and last terminal residues. HYDRO entries are in a 7 to 8 column
 * format: 1) name of hydrogen atom, 2) number of hydrogens needed, 3) bonding
 * type (see add_hydrogens in hbuild.c), 4) distance heavy-hydrogen atom,
 * 5-8) reference atoms for position calculations.  HEAVY entries list the
 * heavy atoms of a residue.  The residue record is terminated with END.
 *
 *
 * $Id: top.c 165 2012-06-29 14:41:27Z hhl $
 *
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>

#include "common.h"
#include "top.h"
#include "util/hashtab.h"
#include "util/hashfuncs.h"
#include "util/util.h"


#define TOP_DELIMITER " \t\n"
#define TOP_LINE_LEN 132



/*
 * tospac: change non alpha-numeric characters of PDB atoms to spaces
 *
 * in:  PDB atom string
 *
 */

static void tospac(char *str)
{
  unsigned int i;


  for (i = 0; i < PDB_ATOM_NAME_LEN-1; i++) {
    if (!isalnum((unsigned char) str[i]) ) {
      str[i] = ' ';
    }
  }

  str[i] = '\0';	
}


/*
 * top_read: read a topology database file and convert to internal structure
 *
 * in:  topol hash table to be filled, name of the top file
 * out: topol hash table filled with data from the file
 *
 */

topol_hash *top_read(topol_hash* top_hash, const char *filename)
{
  unsigned int line_cnt = 0, nrec = 0, nheavy = 0;
  unsigned int in_res = 0;
  unsigned int nname = 0, nent = 0, nfields;
  unsigned int table_size = 0;

  bool first;

  char res_type = '@', mol_type = ' ', buffer[TOP_LINE_LEN];
  char first_term[PDB_ATOM_NAME_LEN], last_term[PDB_ATOM_NAME_LEN];
  char atoms[5][PDB_ATOM_NAME_LEN];
  char *bufp, *resn;
  char *heavy_atom, **heavy_atoms = NULL;

  Hashtable *res_table;
  Hashnode *curr_node;

  topol *top, *top_entry;
  topol_hydro *hydrogen, **hydrogens = NULL;

  struct {
    char resName[PDB_RES_NAME_LEN];
    char first[PDB_ATOM_NAME_LEN];
    char last[PDB_ATOM_NAME_LEN];
  } *term_map, *term;

  FILE *topol_stream;



  if (!(topol_stream = fopen(filename, "r")) ) {
    perror(filename);
    exit(EXIT_FAILURE);
  }

  top = NULL;
  term_map = NULL;

  while (fgets(buffer, TOP_LINE_LEN, topol_stream) ) {  /* read lines */
    line_cnt++;

    if (!(bufp = normln(buffer)) )
      continue;

    if (*bufp == '[') {
      bufp++;

      if (STRNEQ(bufp, "proteins", 8) ) {
	mol_type = 'P';
      } else if (STRNEQ(bufp, "DNA", 3) ) {
	mol_type = 'D';
      } else if (STRNEQ(bufp, "RNA", 3) ) {
	mol_type = 'R';
      } else if (STRNEQ(bufp, "other", 5) ) {
	mol_type = 'O';
      } else {
	prerror(1, "%s: unknown molecule type %s in line %d.\n", filename,
		bufp, line_cnt);
      }
    } else if (STRNEQ(bufp, "RESIDUE", 7) ) {
      if (mol_type == ' ') {
	prerror(1, "%s: no molecule type set in line %d.\n", filename,
		line_cnt);
      }

      if (in_res) {
	prerror(1, "%s: previous residues no properly ENDed in line %d.\n",
		filename, line_cnt);
      }

      bufp += 7;

      nent = 0;
      hydrogens = NULL;

      nheavy = 0;
      heavy_atoms = NULL;

      *first_term = '\0';
      *last_term = '\0';

      nname = 0;

      if (!(bufp = strtok(bufp, TOP_DELIMITER) ) ) {
	prerror(1, "%s: invalid residue entry in line %d.\n", filename, line_cnt);
      }

      while (bufp) {
	nrec++;

	top = reallocate(top, (nrec+1) * sizeof(*top) );

	// FIXME: the PGI optimizer turns top[nrec-1].resName into an invalid
	// address, but why?!
	if (!pdb_format_residue(top[nrec-1].resName, bufp) )
	  prerror (2, "%s: atom name %s too long.\n", filename, bufp);

	top[nrec-1].res_type = res_type;
	top[nrec-1].mol_type = mol_type;
	top[nrec-1].first_term = NULL;
	top[nrec-1].last_term = NULL;

	term_map = reallocate(term_map, (nrec+1) * sizeof(*term_map) );

	if (!pdb_format_residue(term_map[nrec-1].resName, bufp) )
	  prerror (2, "%s: atom name %s too long.\n", filename, bufp);

	term_map[nrec-1].first[0] = '\0';
	term_map[nrec-1].last[0] = '\0';

	bufp = strtok(NULL, TOP_DELIMITER);
	nname++;
      }

      in_res = 1;

      continue;
    } else if (STRNEQ(bufp, "RTYPE", 5) ) {
      bufp += 5;

      nfields = sscanf(bufp, " %c", &res_type);

      if (res_type != 'A' && res_type != 'H' && res_type != '@')
	prerror(2, "%s: unkown residue type found in line %d.\n",
		filename, nfields, line_cnt);

    } else if (STRNEQ(bufp, "FTERM", 5) || STRNEQ(bufp, "LTERM", 5) ) {
      if (!in_res) {
	prerror(1, "%s: not inside residue entry in line %d.\n", filename, line_cnt);
      }

      first = false;

      if (*bufp == 'F') {
	first = true;
      }

      bufp += 5;

      if (!(bufp = strtok(bufp, TOP_DELIMITER) ) ) {
	prerror(1, "%s: invalid first terminal entry in line %d.\n", filename, line_cnt);
      }

      if (first) {
	resn = pdb_format_residue(first_term, bufp);
      } else {
	resn = pdb_format_residue(last_term, bufp);
      }

      if (!resn)
	prerror(1, "%s: residue name cannot be longer than %d characters (line %d).\n",
		filename, PDB_RES_NAME_LEN-1, line_cnt);
    } else if (STRNEQ(bufp, "HYDRO", 5) )  {
      if (mol_type == ' ') {
	prerror(1, "%s: no molecule type set in line %d.\n", filename, line_cnt);
      } else if (!in_res) {
	prerror(1, "%s: not inside residue record in line %d.\n", filename, line_cnt);
      }

      bufp += 5;
      nent++;

      hydrogens = reallocate(hydrogens, (nent+1) * sizeof(*hydrogens));
      hydrogen = allocate(sizeof(*hydrogen));

      *atoms[4] = '\0';

      nfields = sscanf(bufp, "%s %u %u %f %s %s %s %s", atoms[0],
		       &(hydrogen->nhyd), &(hydrogen->type), &(hydrogen->xhdist),
		       atoms[1], atoms[2], atoms[3], atoms[4]);

      if (nfields < 5)
	prerror(2, "%s: only %d fields read successfully in line %d.\n",
		filename, nfields, line_cnt);

      /* format atom names as in PDB entries */
      for (unsigned int i = 0; i < 5; i++) {
	if (*atoms[i] == '<') {
	  strncpy(hydrogen->atoms[i], atoms[i]+1, PDB_ATOM_NAME_LEN-1);
	  hydrogen->atoms[i][PDB_ATOM_NAME_LEN-1] = '\0';

	  tospac(hydrogen->atoms[i]);
	} else {
	  if (pdb_format_atom(hydrogen->atoms[i], atoms[i]) )
	    prerror(2, "%s: atom name %s too long.\n", filename, atoms[i]);
	}
      }

      hydrogens[nent-1] = hydrogen;
    } else if (STRNEQ(bufp, "HEAVY", 5) )  {
      if (!in_res) {
	prerror(1, "%s: not inside residue record in line %d.\n", filename, line_cnt);
      }

      bufp += 5;

      if (!(bufp = strtok(bufp, TOP_DELIMITER) ) ) {
	prerror(1, "%s: invalid heavy atom record in line %d.\n", filename, line_cnt);
      }

      while (bufp) {
	nheavy++;

	heavy_atoms = reallocate(heavy_atoms, (nheavy+1) * sizeof(*heavy_atoms) );
	heavy_atom = allocate(PDB_RES_NAME_LEN);

	if (*bufp == '<') {
	  strncpy(heavy_atom, bufp + 1, PDB_ATOM_NAME_LEN-1);
	  heavy_atom[PDB_ATOM_NAME_LEN-1] = '\0';

	  tospac(heavy_atom);
	} else {
	  if (pdb_format_atom(heavy_atom, bufp) ) {
	    prerror (2, "%s: atom name %s too long.\n", filename, bufp);
	  }
	}

	heavy_atoms[nheavy-1] = heavy_atom;

	bufp = strtok(NULL, TOP_DELIMITER);
      }
    } else if (STRNEQ(bufp, "END", 3) )  {
      if (!in_res) {
	prerror(1, "%s: not inside residue record in line %d.\n", filename, line_cnt);
      }

      if (!hydrogens) {
	prerror(1, "%s: no hydrogen entries found in residue %s (line %d).\n",
		filename, top[nrec-1].resName, line_cnt);
      }

      if (!heavy_atoms) {
	prerror(1, "%s: no heavy atom entries found in residue %s (line %d).\n",
		filename, top[nrec-1].resName, line_cnt);
      }

      top[nrec-1].res_type = res_type;
      heavy_atoms[nheavy] = NULL;
      hydrogens[nent] = NULL;

      for (unsigned int i = nrec - nname; i < nrec; i++) {
	strncpy(term_map[i].first, first_term, PDB_ATOM_NAME_LEN-1);
	term_map[i].first[PDB_ATOM_NAME_LEN-1] = '\0';

	strncpy(term_map[i].last, last_term, PDB_ATOM_NAME_LEN-1);
	term_map[i].last[PDB_ATOM_NAME_LEN-1] = '\0';

	top[i].heavy_atoms = heavy_atoms;
	top[i].hydrogens = hydrogens;
      }

      in_res = 0;
    } else {
      prerror(1, "%s: unknown keyword in line %d.\n", filename, line_cnt);
    }
  }

  fclose(topol_stream);

  if (in_res)
    prerror(1, "%s: last END missing.\n", filename);

  top[nrec].mol_type = '\0';
  term_map[nrec].resName[0] = '\0';

  /*
   * hash residue names and link FTERM and LTERM to their respective entries
   */

  table_size = hibit(nrec + 1) << 1;

  if ( !(res_table = hash_init(&kandr2_hash, table_size) ) ) {
    prerror(1, "Cannot allocate enough memory for hash table.\n");
  }

  for (unsigned int i = 0; i < nrec; i++) {
    hash_insert(res_table, top[i].resName, strlen(top[i].resName), top + i);
  }

  for (term = term_map; *term->resName; term++) {
    if ( !(curr_node = hash_search(res_table, term->resName,
				   strlen(term->resName) ) ) ) {
      prerror(1, "Fatal: key %s not found in hash table.\n", term->resName);
    }

    top_entry = hash_node_get_data(curr_node);

    if (*term->first) {
      if ( !(curr_node = hash_search(res_table, term->first,
				     strlen(term->first) ) ) ) {
	prwarn("%s: terminal residue %s does not exist in topology "
	       "database.\n", filename, term->first);
	top_entry->first_term = NULL;
      } else {
	top_entry->first_term = hash_node_get_data(curr_node);
      }
    }
      
    if (*term->last) {
      if ( !(curr_node = hash_search(res_table, term->last,
				     strlen(term->last) ) ) ) {
	prwarn("%s: terminal residue %s does not exist in topology "
	       "database.\n", filename, term->last);
	top_entry->last_term = NULL;
      } else {
	top_entry->last_term = hash_node_get_data(curr_node);
      }
    }
  }

  free(term_map);
  term_map = NULL;

  top_hash = allocate(sizeof(*top_hash) );
  top_hash->data = top;
  top_hash->hash_table = res_table;

  return top_hash;
}


/*
 * top_print: print a topology database to stdout (simpler format, for debugging)
 *
 * in:  top structure
 *
 */

#ifndef NDEBUG

void top_print(topol* top)
{
  char **h;
  topol *p;
  topol_hydro **es;


  for (p = top; p->mol_type; p++) {
    printf("RESIDUE %s\n", p->resName);

    if (p->first_term)
      printf("  FTERM %s\n", p->first_term->resName);

    if (p->last_term)
      printf("  LTERM %s\n", p->last_term->resName);

    for (es = p->hydrogens; *es; es++) {
      printf("  HYDRO %-4s  %i  %i  %5.3f  %-4s  %-4s  %-4s  %-4s\n", (*es)->atoms[0],
	     (*es)->nhyd, (*es)->type, (*es)->xhdist, (*es)->atoms[1],
	     (*es)->atoms[2], (*es)->atoms[3],
	     (*es)->atoms[4][0] != '\0' ? (*es)->atoms[4] : "");
    }

    printf("  HEAVY");

    for (h = p->heavy_atoms; *h; h++) {
      printf(" %s", *h);
    }

    printf("\n"
	   "END\n\n");
  }
}

#endif


/*
 * top_destroy: release memory allocated for the internal topology database
 *
 * in:  top structure
 *
 */

void top_destroy(topol_hash *top)
{
  unsigned int j, n;

  bool found;

  topol *p;
  topol_hydro **es, ***m;


  for (n = 0, p = top->data; p->mol_type; n++, p++);

  m = allocate(n * sizeof(***m) );
  memset(m, 0, n * sizeof(***m) );
  
  for (j = 0, p = top->data; p->mol_type; p++) {

    found = false;

    for (unsigned int i = 0; i < j; i++) {
      if (p->hydrogens == m[i]) {
	found = true;
	break;
      }
    }

    if (!found) {
      for (es = p->hydrogens; *es; es++) {
	free(*es);
      }

      free(p->hydrogens);
      m[j] = p->hydrogens;
      j++;

      for (char **h = p->heavy_atoms; *h; h++) {
	free(*h);
      }

      free(p->heavy_atoms);
    }
  }

  free(top->data);
  free(m);

  hash_destroy(top->hash_table);

  free(top);
}

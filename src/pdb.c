/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Read and write a file in Protein Data Bank (PDB) format (for standard
 * documents see http://www.wwpdb.org/docs.html).  Also some helper functions
 * are provided here.
 *
 * ATOM/HETATM records must follow the standard closely although a slightly
 * more relaxed regime is used.  SSBOND records are read if desired, one chosen
 * MODEL can be read, TER and CRYST1 records are read.  Some of the title
 * section records and some of the REMARKs are read for output only.
 *
 *
 * $Id: pdb.c 162 2012-06-25 14:33:29Z hhl $
 *
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <float.h>

#include "common.h"
#include "pdb.h"
#include "util/stack.h"
#include "util/util.h"
#include "util/zio.h"


#define PDB_REC_LEN 6
#define PDB_STD_IN_FORMAT "%*6c%5c%*c%4c%c%4c%c%4i%c%*3c%8f%8f%8f%6f%6f%*6c%4c%2c%2c"
#define PDB_STD_OUT_FORMAT "%6s%5s %4s%c%4s%c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n"
#define PDB_MIN_OUT_FORMAT "%6s%5s %4s %4s%c%4i    %8.3f%8.3f%8.3f\n"
#define PDB_SSBOND_FORMAT "%*7c%3i%*5c%c%*c%3i%c%*7c%c%*c%3i%c%*23c%6c%*c%6c%*2c%5f"


enum pdb_format_t {PDB_FMT_STD, PDB_FMT_MIN};



/*
 * pdb_add_atom_node: add a atom node to a link list
 *
 * in:  previous node
 * out: new node
 *
 */

static pdb_atom *pdb_add_atom_node (pdb_atom *old)
{
  pdb_atom *new;


  new = allocate(sizeof(*new));
  memset(new, 0, sizeof(*new));

  new->next = NULL;

  if (old) old->next = new;

  return new;
}


/*
 * pdb_insert_atom_node: insert an atom node between two nodes
 *
 * in:  previous node
 * out: new node
 */

pdb_atom *pdb_insert_atom_node (pdb_atom *old)
{
  pdb_atom *new;


  new = pdb_add_atom_node(NULL);

  new->next = old->next;
  old->next = new;

  return new;
}


/*
 * pdb_add_residue_node: add a residue node to a link list
 *
 * in:  previous node
 * out: new node
 *
 */

static pdb_residue *pdb_add_residue_node (pdb_residue *old)
{
  pdb_residue *new;


  new = allocate(sizeof(*new));
  memset(new, 0, sizeof(*new));

  new->next = NULL;

  if (old) old->next = new;

  return new;
}


/*
 * pdb_add_chain_node: add a chain node to a link list
 *
 * in:  previous node
 * out: new node
 *
 */

static pdb_chain *pdb_add_chain_node (pdb_chain *old)
{
  pdb_chain *new;


  new = allocate(sizeof(*new));
  memset(new, 0, sizeof(*new));

  new->next = NULL;

  if (old) old->next = new;

  return new;
}


/*
 * pdb_format_atom: format atom entry required for PDB
 *
 * in:  destination string, source string
 * out: 1 if error, 0 if no error
 *
 */

int pdb_format_atom (char *restrict dest, const char *restrict src)
{
  size_t len;


  strcpy(dest, "    ");
  len = strlen(src);		// FIXME: source string may not be \0 terminated

  switch (len) {
  case 0:
    *dest = '\0';
    break;

  case 1:
  case 2:
  case 3:
    strncpy(dest + 1, src, len);
    break;

  case 4:
    strncpy(dest, src, len);
    break;

  default:
    return 1;
  }

  dest[PDB_ATOM_NAME_LEN-1] = '\0';

  return 0;
}


/*
 * pdb_format_resiude: format residue entry required for PDB
 *
 * in:  destination string, source string
 * out: destination string or NULL if error
 *
 */

char *pdb_format_residue(char *restrict dest, const char *restrict src)
{
  size_t len;


  len = strlen(src);		// FIXME: source string may not be \0 terminated

  switch (len) {
  case 0:
    return NULL;
    break;

  case 1:
    dest[2] = src[0];
    dest[0] = dest[1] = dest[3] = ' ';
    break;
	
  case 2:
    dest[2] = src[1];
    dest[1] = src[0];
    dest[0] = dest[3] = ' ';
    break;

  case 3:
    strncpy(dest, src, len);
    dest[3] = ' ';
    break;

  case 4:
    strncpy(dest, src, len);
    break;

  default:    
    return NULL;
  }

  dest[PDB_RES_NAME_LEN-1] = '\0';

  return dest;
}


/*
 * extrdat: move past white-space of a string, remove trailing white-space
 *
 * in:  string, offset
 * out: pointer to first non-space
 *
 */

static char *extrdat(char *pos)
{
  char *end;

  pos = pos + strspn(pos, " ");
  end = pos + strlen(pos);

  while (isspace((unsigned char) *--end) && end >= pos) {
    *end = '\0';
  }

  return pos;
}



/*
 * pdb_read: read and analyse ATOM/HETATM, SSBOND, TER, and CRYST1 records from
 *           a file
 *
 * in:  pdb structure, PDB file name, name for CYS residues in disulfide bond,
 *      alternate location indicator, chosen model number, counter for S-S bonds
 * out: pdb root structure
 *
 */

pdb_root *pdb_read(pdb_root *pdb, const char *filename, const char* ss_name,
		   int model_no, int *nssb)
{
  int atom_cnt = 0, residue_cnt = 0, chain_cnt = 0, line_cnt = 0;
  int resSeq, old_resSeq = INT_MIN, curr_model_no, gap;
  int serNum, seqNum1, seqNum2, nss = 0, nfields;

  bool new_chain = false, new_residue = false;
  bool ter_found = false, model_found = false, ss_found = false;
  bool mdltyp_found = false, caveat_found = false;
  bool fr465 = false, fr470 = false, fr475 = false, fr480 = false;

  char altLoc, iCode, old_iCode = '\0', chainID;
  char old_chainID = '\0', ter_chainID = '\0', curr_rectype = '\0';
  char chainID1, chainID2, icode1, icode2;
  char *pos, *res_name, *end;

  float x, y, z, occupancy, tempFactor;
  float Length, f;

  char serial[PDB_SERIAL_LEN], name[PDB_ATOM_NAME_LEN], segID[PDB_SEG_NAME_LEN];
  char element[PDB_ELEMENT_LEN], charge[PDB_CHARGE_LEN];
  char SymOP1[PDB_SSBOND_SYMOP_LEN], SymOP2[PDB_SSBOND_SYMOP_LEN];
  char resName[PDB_RES_NAME_LEN];
  char buffer[PDB_LINE_LEN];

  FILE *pdb_stream;

  Stack *occ_warn = NULL;

  pdb_atom *old_atom, *curr_atom = NULL;
  pdb_residue *old_residue, *curr_residue = NULL;
  pdb_chain *old_chain, *curr_chain = NULL;

  pdb_ssbond *ssbond, **ssbonds;



  if (!(pdb_stream = fzopen(filename, "r")) ) {
    perror(filename);
    exit(2);
  }

  pdb = allocate(sizeof(*pdb) );
  pdb->model_no = 0;
  pdb->ID[0] = '\0';
  pdb->cryst1[0] = '\0';
  pdb->ssbonds = NULL;
  pdb->first_chain = NULL;

  *nssb = 0;

  occ_warn = stack_init(occ_warn);


  while (fzgets(pdb_stream, buffer, PDB_LINE_LEN) ) {
    line_cnt++;
    delnl(buffer);		/* fgets stores the newline */

    if (STRNEQ(buffer, "ATOM", 4) || STRNEQ(buffer, "HETATM", 6)) {
      if (model_found && curr_model_no != model_no) {
	continue;
      }

      serial[0] = name[0] = resName[0] = segID[0] = element[0] =
	charge[0] = '\0';
      altLoc = chainID = iCode = ' ';
      resSeq = 0;
      x = y = z = occupancy = tempFactor = 0.0;

      nfields = sscanf(buffer, PDB_STD_IN_FORMAT,
		       serial, name, &altLoc, resName, &chainID, &resSeq,
		       &iCode, &x, &y, &z, &occupancy, &tempFactor,
		       segID, element, charge);

      if (nfields < 10) {
	prerror(2, "%s: only %d ATOM/HETATM fields read successfully in line %d"
		"but expected at least 10.\n", filename, nfields, line_cnt);
      }

      if (ishydrogen(element, name) ) {
	if (options.remh) {
	  continue;
	} else {
	  strcpy(element, " H"); // tag hydrogens
	}
      }

      resName[PDB_RES_NAME_LEN-1] = '\0';

      if (chainID != old_chainID ||
	  (ter_found && ter_chainID != old_chainID) ) {
	old_chain = curr_chain;
	curr_chain = pdb_add_chain_node(old_chain);

	curr_chain->chainID = chainID;

	if (atom_cnt == 0) {
	  pdb->first_chain = curr_chain;
	}

	new_chain = true;
      }

      if (new_chain || resSeq != old_resSeq ||
	  iCode != old_iCode) {	// if new res
	if (options.warnocc && !stack_is_empty(occ_warn) ) {
	  prwarn("very low occupancy for atoms in residue %s %d%c %c: ",
		 curr_residue->resName, curr_residue->resSeq,
		 curr_residue->iCode, curr_chain->chainID);

	  while ( (res_name = stack_pop(occ_warn)) ) {
	    fprintf(stdout, " %s", (char *)res_name);
	  }

	  fprintf(stdout, "\n");
	}

	old_residue = curr_residue;

	curr_residue = pdb_add_residue_node(old_residue);

	curr_residue->resSeq = resSeq;
	curr_residue->iCode = iCode;

	curr_residue->rectype = curr_rectype = buffer[0];

	strncpy(curr_residue->segID, segID, PDB_SEG_NAME_LEN-1);
	curr_residue->segID[PDB_SEG_NAME_LEN-1] = '\0';

	curr_residue->chain = curr_chain;

	gap = resSeq - old_resSeq - 1;

	if (gap > 0 && !new_chain && buffer[0] == 'A') {
	  prwarn("gap of %i residue%s prior to %s %d%c %c\n",
		 gap, gap > 1 ? "s" : "", resName,  curr_residue->resSeq,
		 curr_residue->iCode, curr_chain->chainID);
	}

	if (STRNEQ(resName, "CYS ", PDB_ATOM_NAME_LEN-1)) {
	  (*nssb)++;

	  if (options.rssb && pdb->ssbonds) {
	    for (ssbonds = pdb->ssbonds; *ssbonds; ssbonds++) {
	      ssbond = *ssbonds;

	      ss_found = false;

	      if ( (resSeq == ssbond->ss1.seqNum &&
		    chainID == ssbond->ss1.chainID) ||
		   (resSeq == ssbond->ss2.seqNum  &&
		    chainID == ssbond->ss2.chainID) ) {
		ss_found = true;
		break;
	      }
	    }

	    if (ss_found) {
	      strncpy(resName, ss_name, PDB_RES_NAME_LEN-1);
	      resName[PDB_RES_NAME_LEN-1] = '\0';
	    }
	  }
	}

	strncpy(curr_residue->resName, resName, PDB_RES_NAME_LEN-1);
	curr_residue->resName[PDB_RES_NAME_LEN-1] = '\0';

	new_residue = true;
      }

      if (!STRNEQ(curr_residue->resName, resName, PDB_RES_NAME_LEN-1) &&
	  iCode == curr_residue->iCode) {
	prerror(1, "residue %s %d%c %c has also other name: %s, "
		"check SEQADV/REMARK 999.\n",
		curr_residue->resName, curr_residue->resSeq, curr_residue->iCode,
		curr_chain->chainID, resName);
      }

      if (buffer[0] != curr_rectype) {
	prerror(1, "residue %s %d%c %c has both ATOM and HETATM records.\n",
		curr_residue->resName, curr_residue->resSeq, curr_residue->iCode,
		curr_chain->chainID);
      }

      old_atom = curr_atom;
      curr_atom = pdb_add_atom_node(old_atom);

      curr_atom->altLoc = altLoc;
      vecCreate(curr_atom->pos, x, y, z);
      curr_atom->occupancy = occupancy;
      curr_atom->tempFactor = tempFactor;

      strncpy(curr_atom->serial, serial, PDB_SERIAL_LEN-1);
      curr_atom->serial[PDB_SERIAL_LEN-1] = '\0';

      strncpy(curr_atom->name, name, PDB_ATOM_NAME_LEN-1);
      curr_atom->name[PDB_ATOM_NAME_LEN-1] = '\0';

      strncpy(curr_atom->element, element, PDB_ELEMENT_LEN-1);
      curr_atom->element[PDB_ELEMENT_LEN-1] = '\0';

      strncpy(curr_atom->charge, charge, PDB_CHARGE_LEN-1);
      curr_atom->charge[PDB_CHARGE_LEN-1] = '\0';

      curr_atom->residue = curr_residue;

      if (options.warnocc && occupancy < FLT_EPSILON) {
	stack_push_uniq(occ_warn, curr_atom->name, PDB_ATOM_NAME_LEN-1);
      }

      /* FIXME: check for "broken" chains! */
      if (new_chain) {
	curr_chain->first_residue = curr_residue;

	chain_cnt++;
	old_chainID = chainID;

	ter_found = false;
	new_chain = false;
      }

      if (new_residue) {
	curr_residue->first_atom = curr_atom;

	residue_cnt++;
	old_resSeq = resSeq;
	old_iCode = iCode;

	new_residue = false;
      }

      atom_cnt++;
    } else if (STRNEQ(buffer, "MODEL", 5) ) {
      model_found = true;
	
      nfields = sscanf(buffer, "%*10c%4i", &curr_model_no);

      if (nfields != 1) {
	prerror(2, "%s: MODEL record requires serial in line %d.\n",
		filename, line_cnt);
      }

      if (model_no < 0) {
	model_no = curr_model_no;
      }
    } else if (STRNEQ(buffer, "TER", 3) ) {
      ter_found = true;
      nfields = sscanf(buffer, "%*21c%c", &ter_chainID);

      if (nfields != 1) {
	ter_chainID = ' ';
      }
    } else if (options.rssb && STRNEQ(buffer, "SSBOND", 6) ) {
      SymOP1[0] = SymOP2[0] = '\0';
      Length = 0.0;

      nfields = sscanf(buffer, PDB_SSBOND_FORMAT,
		       &serNum, &chainID1, &seqNum1, &icode1,
		       &chainID2, &seqNum2, &icode2, SymOP1, SymOP2, &Length);

      if (nfields < 9) {
	prwarn("%s: only %d SSBOND fields read successfully in line %d but"
	       "expected at least 10.\n", filename, nfields, line_cnt);
      }

      nss++;
      pdb->ssbonds = reallocate(pdb->ssbonds, (nss+1) * sizeof(*pdb->ssbonds));
      ssbond = allocate(sizeof(*ssbond));

      ssbond->serNum = serNum;

      ssbond->ss1.chainID = chainID1;
      ssbond->ss1.seqNum = seqNum1;
      ssbond->ss1.icode = icode1;

      ssbond->ss2.chainID = chainID2;
      ssbond->ss2.seqNum = seqNum2;
      ssbond->ss2.icode = icode2;

      strncpy(ssbond->ss1.SymOP, SymOP1, PDB_SSBOND_SYMOP_LEN-1);
      ssbond->ss1.SymOP[PDB_SSBOND_SYMOP_LEN-1] = '\0';

      strncpy(ssbond->ss2.SymOP, SymOP2, PDB_SSBOND_SYMOP_LEN-1);
      ssbond->ss2.SymOP[PDB_SSBOND_SYMOP_LEN-1] = '\0';

      ssbond->Length = Length;

      pdb->ssbonds[nss-1] = ssbond;
      pdb->ssbonds[nss] = NULL;
    } else if (STRNEQ(buffer, "CRYST1", 6) ) {
      strncpy(pdb->cryst1, buffer, PDB_LINE_LEN-1);
      pdb->cryst1[PDB_LINE_LEN-1] = '\0';
    } else if (STRNEQ(buffer, "HEADER", 6) )  {
      prnote("header and title of %s\n   %s\n", filename, extrdat(buffer + 6) );

      if (strlen(buffer) > 66) {
	strncpy(pdb->ID, buffer + 62, PDB_ID_LEN-1);
	pdb->ID[PDB_ID_LEN-1] = '\0';
      }
    } else if (STRNEQ(buffer, "OBSLTE", 6) )  {
      prwarn("this PDB has been obsoleted by %s\n", extrdat(buffer+31) );
    } else if (STRNEQ(buffer, "TITLE", 5) )  {
      fprintf(stdout, "   %s\n", extrdat(buffer+10) );
    } else if (STRNEQ(buffer, "SPLIT", 5) )  {
      prwarn("PDB has been split.  Required IDs to reconstitute: %s\n",
	     extrdat(buffer+11) );
    } else if (STRNEQ(buffer, "CAVEAT", 6) )  {
      if (!caveat_found) {
	prwarn("This PDB contains SEVERE ERRORS:\n");
	caveat_found = true;
      }

      fprintf(stdout, "    %s\n", extrdat(buffer+10) );
    } else if (STRNEQ(buffer, "EXPDTA", 6) )  {
      prnote("PDB reports experiment type as %s\n", extrdat(buffer+6) );
    } else if (STRNEQ(buffer, "NUMMDL", 6) )  {
      pos = buffer + 10;
      pos[14] = '\0';
      prnote("PDB contains %i models\n", atoi(pos));
    } else if (STRNEQ(buffer, "MDLTYP", 5) )  {

      if (!mdltyp_found) {
	prnote("PDB reports model type as\n");
	mdltyp_found = true;
      }

      fprintf(stdout, "    %s\n", extrdat(buffer+10) );
    } else if (STRNEQ(buffer, "REMARK   2 RESOLUTION.", 22)) {
      pos = buffer + 22;
      errno = 0;
      f = strtof(pos, &end);

      if ( !(end == pos || errno == ERANGE) ) {
	prnote("PDB resolution is %.2f\n", f);
      }
    } else if (STRNEQ(buffer, "REMARK   4 ", 11) ) {
      if (STRNEQ(buffer + 30, "FORMAT", 6)) {
	prnote("PDB version %s\n", extrdat(buffer+40));
      }
    } else if (STRNEQ(buffer, "REMARK 2", 8) && 
	       (buffer[8] == '0' || buffer[8] == '1' || buffer[8] == '3' ||
		buffer[8] == '4' || buffer[8] == '5' || buffer[8] == '6') &&
	       STRNEQ(buffer + 11, " PH ", 4) ) {
      pos = strchr(buffer, ':');

      if (pos) {
	pos++;
	errno = 0;
	f = strtof(pos, &end);

	if ( !(end == pos || errno == ERANGE) ) {
	  prnote("PDB reports a pH of %.2f in REMARK 2nn\n", f);
	}
      }
    } else if (STRNEQ(buffer, "REMARK 465", 10) ) {
      if (!fr465) {
	prnote("PDB warns of missing residues\n");
	fr465 = true;
      }
    } else if (STRNEQ(buffer, "REMARK 470", 10) ) {
      if (!fr470) {
	prnote("PDB warns of missing atoms\n");
	fr470 = true;
      }
    } else if (STRNEQ(buffer, "REMARK 475", 10) ) {
      if (!fr475) {
	prnote("PDB warns of residues with zero occupancy\n");
	fr475 = true;
      }
    } else if (STRNEQ(buffer, "REMARK 480", 10) ) {
      if (!fr480) {
	prnote("PDB warns of non-hydrogens with zero occupancy\n");
	fr480 = true;
      }
    } 
  }

  fzclose(pdb_stream);

  if (!pdb->first_chain) {
    prerror(2, "\n%d lines read but no atoms extracted from %s.\n",
	    line_cnt, filename);
  }

  if (model_found) {
    pdb->model_no = model_no;
  }

  /* flush warning of last residue if it exists */
  if (options.warnocc && !stack_is_empty(occ_warn) ) {
    prwarn("very low occupancy for atoms in residue %s %d%c %c: ",
	   curr_residue->resName, curr_residue->resSeq, curr_residue->iCode,
	   curr_chain->chainID);

    while ( (res_name = stack_pop(occ_warn)) ) {
      fprintf(stdout, " %s", (char *)res_name);
    }

    fprintf(stdout, "\n");
  }

  stack_destroy(occ_warn);

  fprintf(stdout, "\n%d atoms, %d residues, %d chain%s read\n",
	  atom_cnt, residue_cnt, chain_cnt, chain_cnt > 1 ? "s" : "");

  pdb->natoms = atom_cnt;
  pdb->nres = residue_cnt;
  pdb->nchains = chain_cnt;


  return pdb;
}


/*
 * pdb_write: write a PDB file in either standard or relaxed standard format
 *
 * in:  pdb root structure, file name, chosen format, name of CYS residue in
 *      disulfide bond
 *
 */

void pdb_write(pdb_root *pdb, const char *filename, const char *format,
	       const char* ss_name, char altLoc)
{
  int std_type= PDB_FMT_STD, resSeq = 0;
  int serno = 0, atom_cnt = 0, residue_cnt = 0, chain_cnt = 0;

  char chainID = ' ', iCode = ' ';
  char *resName = NULL, *rectype = NULL;
  char serial[6];
  char atomrec[] = "ATOM  ", hetrec[] = "HETATM";

  FILE *pdb_stream;

  pdb_atom *curr_atom;
  pdb_residue *curr_residue;
  pdb_chain *curr_chain;

  pdb_ssbond *ssbond, **ssbonds;


  if (*format == '\0' || STRNEQ(format, "std", 3) ) {
    std_type = PDB_FMT_STD;
    options.noter = 0;
  } else if STRNEQ(format, "min", 3) {
    std_type = PDB_FMT_MIN;
  } else {
    prerror(2, "Unknown format type: %s\n", format);
  }

  if (!(pdb_stream = fopen(filename, "w")) ) {
    perror(filename);
    exit(2);
  }

  if (pdb->ID[0] != '\0') {
    fprintf(pdb_stream, "REMARK   this is a conversion of PDB ID %s\n",
	    pdb->ID);
  }

  if (options.wrss && pdb->ssbonds) {
    for (ssbonds = pdb->ssbonds; *ssbonds; ssbonds++) {
	ssbond = *ssbonds;

      fprintf(pdb_stream,
	      "SSBOND %3i CYS %c %4i%c   CYS %c %4i%c                       "
	      "%6s %6s",
	      ssbond->serNum, ssbond->ss1.chainID, ssbond->ss1.seqNum,
	      ssbond->ss1.icode, ssbond->ss2.chainID, ssbond->ss2.seqNum,
	      ssbond->ss2.icode, ssbond->ss1.SymOP, ssbond->ss2.SymOP);

      if (ssbond->Length > 0.0) {
	fprintf(pdb_stream, " %5.2f\n", ssbond->Length);
      } else {
	fprintf(pdb_stream, "      \n");
      }
    }
  }

  if (pdb->cryst1[0] != '\0' && !options.nocryst) {
    switch (std_type) {
    case PDB_FMT_STD:
      fprintf(pdb_stream, "%-80s\n", pdb->cryst1);
      break;
    case PDB_FMT_MIN:
      fprintf(pdb_stream, "%s\n", pdb->cryst1);
      break;
    }
  }

  if (pdb->model_no > 0 && !options.nomodel) {
    switch (std_type) {
    case PDB_FMT_STD:
      fprintf(pdb_stream, "MODEL     %4i%65s\n", pdb->model_no, " ");
	break;
    case PDB_FMT_MIN:
      fprintf(pdb_stream, "MODEL     %4i\n", pdb->model_no);
      break;
    }
  }

  for (curr_chain = pdb->first_chain; curr_chain; curr_chain = curr_chain->next) {
    chain_cnt++;		/* chain */

    for (curr_residue = curr_chain->first_residue;
	 curr_residue && curr_residue->chain == curr_chain;
	 curr_residue = curr_residue->next) {  /* residue */
      residue_cnt++;

      if (curr_residue->rectype == 'A') {
	rectype = atomrec;
      } else if (curr_residue->rectype == 'H') {
	rectype = hetrec;
      } else {
	prerror(1, "rectype %c is unknown.\n", curr_residue->rectype);
      }


      for (curr_atom = curr_residue->first_atom;
	   curr_atom && curr_atom->residue == curr_residue;
	   curr_atom = curr_atom->next) {  /* atom */

	if (curr_atom->altLoc != altLoc && curr_atom->altLoc != ' ') {
	  continue;
	}

	atom_cnt++;

	if (options.keepser) {
	  strncpy(serial, curr_atom->serial, PDB_SERIAL_LEN-1);
	  serial[PDB_SERIAL_LEN-1] = '\0';
	} else {
	  serno++;

	  /* FIXME: serial number, also: overflow possible */
	  if (serno > 99999) {
	    serno = 1;
	  }

	  sprintf(serial, "%i", serno);
	}

	if (pdb->ssbonds && !options.keepssn &&
	    STRNEQ(curr_residue->resName, ss_name, PDB_ATOM_NAME_LEN-1) ) {
	  strcpy(curr_residue->resName, "CYS ");
	  curr_residue->resName[PDB_ATOM_NAME_LEN-1] = '\0';
	}

	switch (std_type) {
	case PDB_FMT_STD:
	  fprintf(pdb_stream, PDB_STD_OUT_FORMAT, rectype,
		  serial, curr_atom->name, curr_atom->altLoc,
		  curr_residue->resName, curr_chain->chainID,
		  curr_residue->resSeq, curr_residue->iCode,
		  curr_atom->pos[0], curr_atom->pos[1], curr_atom->pos[2],
		  curr_atom->occupancy, curr_atom->tempFactor,
		  curr_residue->segID, curr_atom->element,
		  curr_atom->charge);
	  break;
	case PDB_FMT_MIN:
	  fprintf(pdb_stream, PDB_MIN_OUT_FORMAT, rectype,
		  serial, curr_atom->name,
		  curr_residue->resName, curr_chain->chainID,
		  curr_residue->resSeq,
		  curr_atom->pos[0], curr_atom->pos[1], curr_atom->pos[2]);
	  break;
	}
      }	/* atom */

      resName = curr_residue->resName;
      chainID = curr_chain->chainID;
      iCode = curr_residue->iCode;
      resSeq = curr_residue->resSeq;
    } /* residue */

    if (rectype[0] == 'A' && !options.noter) {
      switch (std_type) {
      case PDB_FMT_STD:
	if (!options.keepser) {
 	  sprintf(serial, "%i", ++serno);
	}

	fprintf(pdb_stream, "TER   %5s      %4s%c%4i%c%53s\n", serial, resName,
		chainID, resSeq, iCode, " ");
	break;
      case PDB_FMT_MIN:
	fprintf(pdb_stream, "TER\n");
	break;
      }
    }
  } /* chain */

  if (pdb->model_no > 0 && !options.nomodel) {
    switch (std_type) {
    case PDB_FMT_STD:
      fprintf(pdb_stream, "%-80s\n", "ENDMDL");
      break;
    case PDB_FMT_MIN:
      fprintf(pdb_stream, "ENDMDL\n");
      break;
    }
  }
	
  if (!options.noend) {
    switch (std_type) {
    case PDB_FMT_STD:
      fprintf(pdb_stream, "%-80s", "END");
      break;
    case PDB_FMT_MIN:
      fprintf(pdb_stream, "END");
      break;
    }
  }

  fprintf(stdout, "%d atoms, %d residues, %d chain%s written\n",
	  atom_cnt, residue_cnt, chain_cnt, chain_cnt > 1 ? "s" : "");

  fclose(pdb_stream);
}


/*
 * pdb_destroy: release allocated memory associated with pdb structures
 *
 * in:  pdb root structure
 *
 */

void pdb_destroy(pdb_root *pdb)
{
  pdb_atom *curr_atom, *next_atom;
  pdb_residue *curr_residue, *next_residue;
  pdb_chain *curr_chain, *next_chain;

  pdb_ssbond **ssbonds;


  if (pdb->ssbonds) {
    for (ssbonds = pdb->ssbonds; *ssbonds; ssbonds++) {
      free(*ssbonds);
    }
  }

  free(pdb->ssbonds);

  for (curr_chain = pdb->first_chain; curr_chain; curr_chain = next_chain) {
    next_chain = curr_chain->next;

    for (curr_residue = curr_chain->first_residue;
	 curr_residue && curr_residue->chain == curr_chain;
	 curr_residue = next_residue) {
      next_residue = curr_residue->next;

      for (curr_atom = curr_residue->first_atom;
	   curr_atom && curr_atom->residue == curr_residue;
	   curr_atom = next_atom) {
	next_atom = curr_atom->next;
	free(curr_atom);
      }

      free(curr_residue);
    }

    free(curr_chain);
  }

  free(pdb);
}

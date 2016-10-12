/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Add hydrogens to heavy atoms as defined in a hydrogen database.  Terminal
 * residues of proteins and nucleotides can be treated specially.
 *
 *
 * $Id: hbuild.c 162 2012-06-25 14:33:29Z hhl $
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>

#include "common.h"
#include "pdb.h"
#include "top.h"
#include "util/hashtab.h"
#include "util/stack.h"
#include "util/queue.h"
#include "util/util.h"


#define MAX_XHDIST 1.5		/* "generous" X-H distance squared */
#define ISHYD(e) ( ( (e)[0] ) == ' ' && ( (e)[1] ) == 'H' )

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* pre-calculated value for calculation of H positions */
#define SIN_tetra 0.9428090415820634       /* sin(109.47) */
#define COS_tetra -0.3333333333333333      /* cos(109.47) */
#define SIN_tetra_h 0.8164965809277260     /* sin(109.47 / 2) */
#define COS_tetra_h 0.5773502691896257     /* cos(109.47 / 2) */
#define SIN_tetra_05 0.4714045207910317	   /* sin(109.47) * 0.5 */
#define SIN_120 0.8660254037844387         /* sin(120) */
#define COS_120 -0.5		           /* cos(120) */

// _very_ arbitrary pre-calculated values for spherical coordinates for H2O
// theta1 = 50 deg, theta2 = 50 + 104.52 deg, phi = 70 deg
#define SIN_theta1 0.76604444311897803520
#define SIN_theta2 0.43019600888661864331
#define COS_phi 0.34202014332566873305
#define SIN_phi 0.93969262078590838405
#define COS_theta1 0.64278760968653932632
#define COS_theta2 0.90273550608028281612



/*
 * fill_atom: fill records of a PDB atom entry
 *
 * in:  atom node, record type, atom name, coordinates
 *
 */

static void fill_atom(pdb_atom *at, const char *atom0, const fvec pos)
{
  strncpy(at->name, atom0, PDB_ATOM_NAME_LEN-1);
  at->name[PDB_ATOM_NAME_LEN-1] = '\0';

  at->altLoc = ' ';
  vecCopy(at->pos, pos);
  at->occupancy = 1.0;
  at->tempFactor = 0.0;

  strcpy(at->element, " H");
  strcpy(at->charge, "  ");
}


/*
 * format_atom_name: append number to PDB atom name
 *
 * in:  name from top (assumed to be PDB_ATOM_NAME_LEN chars long),
 *      character (number) to be added to name
 * out: formatted name
 *
 */

static char *format_atom_name(char *name, char c)
{

  if (name[3] != ' ') {
    for (unsigned int i = 0; i < 3; i++) {
      name[i] = name[i+1];
    }

    name[3] = ' ';
  }

  if (name[2] == ' ') {
    name[2] = c;
  } else if (name[3] == ' ') {
    name[3] = c;
  }

  return name;
}


/*
 * vec_helper: simple helper function for some needed vectors
 *
 * in:  vectors
 *
 */

static void vec_helper(fvec v1, fvec v2, fvec v3,
		       const fvec atom0, const fvec ctrl0, const fvec ctrl1)
{
  vecSub(v1, atom0, ctrl0);
  vecSub(v3, ctrl0, ctrl1);
  vecCrossProd(v2, v1, v3);

  vecScalarDiv(v1, v1, vecLen(v1));
  vecScalarDiv(v2, v2, vecLen(v2));
  vecCrossProd(v3, v2, v1);
}


/*
 * res_check: check if a residue has all heavy atoms as per topology database
 *
 * in:  chain, reside, topology entry
 *
 */

static void res_check(const pdb_chain *chain, const pdb_residue *residue,
		      const topol *entry, char altLoc)
{
  unsigned int n_bb_found = 0, n_CB_found = 0;

  bool found = false;

  char *res_name, **heavy;

  /* FIXME: we should not hardcode this... */
  static const char *const backbone[] = {" CA ", " N  ", " C  ", " O  "};

  Stack *warn = NULL;

  pdb_atom *atom;



  warn = stack_init(warn);

  for (heavy = entry->heavy_atoms; *heavy; heavy++) {
    for (atom = residue->first_atom;
         atom && atom->residue == residue;
         atom = atom->next) {

      if (atom->altLoc != altLoc && atom->altLoc != ' ') {
	continue;
      }

      found = false;

      if (STRNEQ(atom->name, *heavy, PDB_ATOM_NAME_LEN-1) ) {
        found = true;

        /* FIXME: since we do the check for backbone atoms here we could skip
           bb entries in the HEAVY records of the topology database */
        if (entry->mol_type == 'P') {
          for (unsigned int i = 0; i < 4; i++) {
            if (STRNEQ(*heavy, backbone[i],  PDB_ATOM_NAME_LEN-1) ) {
              n_bb_found++;
            }
          }

          if (STRNEQ(*heavy, " CB ",  PDB_ATOM_NAME_LEN-1) ) {
            n_CB_found++;
          }
        }

        break;
      }
    }

    if (!found) {
      stack_push_uniq(warn, *heavy, PDB_RES_NAME_LEN-1);
    }
  }

  if (!stack_is_empty(warn)) {
    prwarn("atoms not found in residue %s %d%c %c: ",
	   residue->resName, residue->resSeq, residue->iCode, chain->chainID);

    while ( (res_name = stack_pop(warn)) ) {
      fprintf(stdout, " %s", (char *)res_name);
    }

    fprintf(stdout, "\n");
  }

  stack_destroy(warn);
}


/*
 * search_top_hydrogens: search for the proper hydrogen entry in the topology
 *
 * in:  hydrogen pointer in topology, atom entry
 * out: pointer to hydrogen entry or NULL if not found
 *
 */

static topol_hydro *search_top_hydrogens(topol_hydro **top_hydro,
					 const char *atom_name)
{
  topol_hydro *entry, **es;


  /* simple linear search */
  for (es = top_hydro; *es; es++) {
    entry = *es;

    if (STRNEQ(entry->atoms[1], atom_name, PDB_ATOM_NAME_LEN-1) ) {
      return entry;
    }
  }

  return NULL;
}


/*
 * search_top_atom: search for the proper control atom in a top entry
 *
 * in:  residue, entry atom name
 * out: entry position
 *
 */

static pdb_atom *search_top_atom(const pdb_residue *residue, const char *name)
{
  pdb_atom *atom;


  for (atom = residue->first_atom;
       atom && atom->residue == residue;
       atom = atom->next) {
    if (ISHYD(atom->element) )
      continue;

    if (STRNEQ(name, atom->name, PDB_ATOM_NAME_LEN-1) ) {
      return atom;
    }
  }

  return NULL;
}


/*
 * add_hydrogens: compute positions for hydrogens according to bonding type
 *
 * in:  atom, current topology entry, current residue, previous residue
 *
 */

static bool add_hydrogens(pdb_atom *atom0,  const topol_hydro *entry,
			  pdb_residue *restrict curr_residue,
			  pdb_residue *restrict prev_residue)
{
  unsigned int ub;

  char name[PDB_ATOM_NAME_LEN];

  float vlen;

  fvec v1, v2, v3, rcent;
  fvec ctrl_atom[3], rH[3];

  pdb_atom *atom, *newH;



  switch (entry->type) {
  case 5:			// 5 contol atoms
    ub = 5;
    break;

  case 10:			// special case water with no control atoms
    ub = 0;
    break;

  default:
    ub = 4;
  }

  for (unsigned int i = 2; i < ub; i++) {
    if (strchr(entry->atoms[i], '-') && prev_residue) {	 // check for prev res
      strncpy(name, entry->atoms[i], PDB_RES_NAME_LEN-1);
      name[PDB_RES_NAME_LEN-1] = '\0';

      if (name[0] == '-') {
	name[0] = ' ';
      } else if (name[1] == '-') {
	name[1] = name[2];
	name[2] = name[3];
	name[3] = ' ';
      }

      atom = search_top_atom(prev_residue, name);

      if (atom) {
	vecCopy(ctrl_atom[i-2], atom->pos);
      } else {
	return false;
      }
    } else {
      atom = search_top_atom(curr_residue, entry->atoms[i]);

      if (atom) {
	vecCopy(ctrl_atom[i-2], atom->pos);
      } else {
	return false;
      }
    }
  }


  switch (entry->type) {
  case 1:			/* planar hydrogens */
    vecSub(v1, atom0->pos, ctrl_atom[0]);
    vecSub(v3, atom0->pos, ctrl_atom[1]);

    vecAdd(v2, v1, v3);
    vecScalarDiv(v2, v2, vecLen(v2));

    for (unsigned int i = 0; i < 3; i++) {
      rH[0][i] = atom0->pos[i] + entry->xhdist * v2[i];
    }

    break;

  case 2:			/* hydrogen bound to O or S */
    vec_helper(v1, v2, v3, atom0->pos, ctrl_atom[0], ctrl_atom[1]);

    for (unsigned int i = 0; i < 3; i++) {
      rH[0][i] = atom0->pos[i] +
	entry->xhdist * SIN_tetra * v3[i] - entry->xhdist * COS_tetra * v1[i];
    }

    break;

  case 3:			/* two planar hydrogens */
    vec_helper(v1, v2, v3, atom0->pos, ctrl_atom[0], ctrl_atom[1]);

    for (unsigned int i = 0; i < 3; i++) {
      rH[0][i] = atom0->pos[i] -
	entry->xhdist * SIN_120 * v3[i] - entry->xhdist * COS_120 * v1[i];
      rH[1][i] = atom0->pos[i] +
	entry->xhdist * SIN_120 * v3[i] - entry->xhdist * COS_120 * v1[i];
    }

    break;

  case 4:			/* three tetrahedal hydrogens */
    vec_helper(v1, v2, v3, atom0->pos, ctrl_atom[0], ctrl_atom[1]);

    for (unsigned int i = 0; i < 3; i++) {
      rH[0][i] = atom0->pos[i] +
	entry->xhdist * SIN_tetra * v3[i]  - entry->xhdist * COS_tetra * v1[i];
      rH[1][i] = atom0->pos[i] -
	entry->xhdist * SIN_tetra_05 * v3[i] +
	entry->xhdist * SIN_tetra_h * v2[i] -
	entry->xhdist * COS_tetra * v1[i];
      rH[2][i] = atom0->pos[i] -
	entry->xhdist * SIN_tetra_05 * v3[i] -
	entry->xhdist * SIN_tetra_h * v2[i] -
	entry->xhdist * COS_tetra * v1[i];
    }

    break;

  case 5:			/* one tetrahedral hydrogen */
    for (unsigned int i = 0; i < 3; i++) {
      rcent[i] = atom0->pos[i] - 
	(ctrl_atom[0][i] + ctrl_atom[1][i] + ctrl_atom[2][i]) / 3.0;
    }

    vlen = vecLen(rcent);

    if (vlen < 0.2) {		// arbitrary length for "short" vector
      for (unsigned int i = 0; i < 3; i++) {
	v1[i] = ctrl_atom[1][i] - ctrl_atom[0][i];
	v2[i] = ctrl_atom[2][i] - ctrl_atom[0][i];
	v3[i] = copysignf(v3[i], rcent[i]);  // gives correct direction?
      }

      vecCrossProd(v3, v1, v2);
      vecScalarDiv(rcent, v3, vecLen(v3));
    } else {
      vecScalarDiv(rcent, rcent, vlen);
    }

    for (unsigned int i = 0; i < 3; i++)
      rH[0][i] = atom0->pos[i] + entry->xhdist * rcent[i];

    break;

  case 6:			/* two tetrahedral hydrogens */
    for (unsigned int i = 0; i < 3; i++) 
      rcent[i] = atom0->pos[i] - (ctrl_atom[0][i] + ctrl_atom[1][i]) / 2.0;

    vecSub(v1, atom0->pos, ctrl_atom[0]);
    vecSub(v2, atom0->pos, ctrl_atom[1]);
    vecCrossProd(v3, v1, v2);

    vecScalarDiv(rcent, rcent, vecLen(rcent));
    vecScalarDiv(v3, v3, vecLen(v3));

    for (unsigned int i = 0; i < 3; i++) {
      rH[0][i] = atom0->pos[i] +
	entry->xhdist * (COS_tetra_h * rcent[i] + SIN_tetra_h * v3[i]);
      rH[1][i] = atom0->pos[i] +
	entry->xhdist * (COS_tetra_h * rcent[i] - SIN_tetra_h * v3[i]);
    }

    break;

    // FIXME: we may want to use a better scheme...
  case 10:			// 3 point water (as TIP3P)
    rH[0][0] = atom0->pos[0] + entry->xhdist * SIN_theta1 * COS_phi;
    rH[0][1] = atom0->pos[1] + entry->xhdist * SIN_theta1 * SIN_phi;
    rH[0][2] = atom0->pos[2] + entry->xhdist * COS_theta1;

    rH[1][0] = atom0->pos[0] + entry->xhdist * SIN_theta2 * COS_phi;
    rH[1][1] = atom0->pos[1] + entry->xhdist * SIN_theta2 * SIN_phi;
    rH[1][2] = atom0->pos[2] - entry->xhdist * COS_theta2;

    break;

  default:
    prerror(1, "hydrogen type %d does not exist in database\n", entry->type);
  }

  for (unsigned int i = 0; i < entry->nhyd; i++) {
    strncpy(name, entry->atoms[0], PDB_RES_NAME_LEN-1);
    name[PDB_RES_NAME_LEN-1] = '\0';

    newH = pdb_insert_atom_node(atom0);
    newH->residue = curr_residue;

    if (entry->nhyd > 1) {
      fill_atom(newH, format_atom_name(name, 48 +  entry->nhyd - i), rH[i]);
    } else {
      fill_atom(newH, name, rH[i]);
    }
  }

  return true;
}


/*
 * hbuild: main loop over heavy atoms and add hydrogens accordingly (actual
 *         routine is in add_hydrogens), check terminal residues and transform
 *         if required
 *
 * in:  pdb root structure, top hash table
 *
 */

void hbuild(pdb_root *pdb, const Hashtable *top, char altLoc)
{
  unsigned int nH;

  bool add_ok;

  char *res_name;

  float dist;

  topol *top_entry;
  topol_hydro *entry;

  Queue *warn = NULL;

  Hashnode *curr_node;

  pdb_atom *atom1, *atom2;
  pdb_residue *curr_residue, *prev_residue;
  pdb_chain *chain;



  warn = queue_init(warn);

  for (chain = pdb->first_chain; chain; chain = chain->next) { /* chain */

    prev_residue = NULL;

    for (curr_residue = chain->first_residue;
	 curr_residue && curr_residue->chain == chain;
	 prev_residue = curr_residue,
	   curr_residue = curr_residue->next) { /* residue */

      if (!(curr_node = hash_search(top, curr_residue->resName,
				    strlen(curr_residue->resName) ) ) ) {
	queue_push_uniq(warn, curr_residue->resName, PDB_RES_NAME_LEN-1);
	continue;
      }

      top_entry = hash_node_get_data(curr_node);

      /* check only according to defined residue type in top.dat as residues
	 may have names colliding with force field conventions, e.g.
	 TYM = TRYPTOPHANYL-5'AMP, HID = (5-HYDROXY-1H-INDOL-3-YL)ACETIC ACID,
	 DGN = D-GLUTAMINE, etc. */
      if (top_entry->res_type != '@' &&
	  curr_residue->rectype != top_entry->res_type) {
	queue_push_uniq(warn, curr_residue->resName, PDB_RES_NAME_LEN-1);
	continue;
      }

      if (curr_residue == chain->first_residue &&
		 top_entry->first_term &&
	  ( (top_entry->mol_type == 'P' && options.nterm) ||
	    (top_entry->mol_type == 'D' && options.dna5term) ||
	    (top_entry->mol_type == 'R' && options.rna5term) ) ) {

	top_entry = top_entry->first_term;
      } else if ( (!curr_residue->next ||
		   curr_residue->next->chain != chain) &&
		  top_entry->last_term &&
		  ( (top_entry->mol_type == 'P' && options.cterm) ||
		    (top_entry->mol_type == 'D' && options.dna3term) ||
		    (top_entry->mol_type == 'R' && options.rna3term) ) ) {

	top_entry = top_entry->last_term;
      }

      res_check(chain, curr_residue, top_entry, altLoc);

      for (atom1 = curr_residue->first_atom;
	   atom1 && atom1->residue == curr_residue;
	   atom1 = atom1->next) { /* atom1 */

	if ( (atom1->altLoc != altLoc && atom1->altLoc != ' ') ||
	     ISHYD(atom1->element) ) {
	  continue;
	}

	nH = 0;

	/* don't look backwards here as we may find badly attached hydrogens,
	   but could make that a check... */
	for (atom2 = atom1->next;
	     atom2 && atom2->residue == curr_residue;
	     atom2 = atom2->next) {

	  if (atom2->altLoc != altLoc && atom2->altLoc != ' ') {
	    continue;
	  }

	  dist = vecDist(atom1->pos, atom2->pos);

	  if (dist < MAX_XHDIST && ISHYD(atom2->element)) {
	    nH++;
	  }
	}

	entry = search_top_hydrogens(top_entry->hydrogens, atom1->name);

	if (entry) {
	  if (nH > entry->nhyd && prev_residue) {
	    prwarn("atom %s-%s %d%c %c has too many hydrogens (%d) already.\n",
		   atom1->name, curr_residue->resName, curr_residue->resSeq,
		   chain->chainID, curr_residue->iCode, nH);
	    continue;
	  } else if (entry->nhyd == nH) {
	    continue;
	  } else if (nH > 0) {
	    if (prev_residue)
	      prwarn("atom %s-%s %d%c %c: cannot handle partially (%d) "
		     "populated hydrogens\n", atom1->name,
		     curr_residue->resName, curr_residue->resSeq,
		     chain->chainID, curr_residue->iCode, nH);

	    continue;
	  } else {
	    add_ok = add_hydrogens(atom1, entry, curr_residue, prev_residue);

	    if (!add_ok) {
	      prwarn("cannot find all control atoms for atom %s (%s %d%c %c) "
		     "in PDB.\n", atom1->name, curr_residue->resName,
		     curr_residue->resSeq, curr_residue->iCode, chain->chainID);
	    }
	  }
	}
      }
    }
  }

  if (!queue_is_empty(warn)) {
    prwarn("residues not found in topology database:");

    while ( (res_name = queue_pop_front(warn)) ) {
      fprintf(stdout, " %s", (char *)res_name);
    }

    fprintf(stdout, "\n");
  }

  queue_destroy(warn);
}

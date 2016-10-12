/*
 * Copyright (C) 2011 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Find S-S bonds and tag "CYS" residues with a user-supplied name if a second
 * "CYS" is within a certain distance.  The code handles only standard "CYS"
 * residues in "ATOM  " records.
 *
 *
 * $Id: ssbuild.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <string.h>
#include <math.h>

#include "common.h"
#include "pdb.h"
#include "util/util.h"


#define MAX_SSDIST 9.0		/*  S-S bond distance */
#define STD_CYS_NAME "CYS "

/* only accept CYS in ATOM records and assume sulfur's first letter is 'S' */
#define NOT_CYS(res,at)							\
  (res->rectype != 'A' ||						\
   !STRNEQ(res->resName, STD_CYS_NAME, PDB_RES_NAME_LEN-1) ||		\
   STRNEQ(res->resName, ss_name, PDB_RES_NAME_LEN-1) ||			\
   !(at->name[1] == 'S' && at->name[2] == 'G') )



/*
 * ssbuild: rename "CYS " residues when in disulfide bond
 *
 * in:  pdb structure, new name for the "CYS " residue
 * out: modified pdb structure
 *
 */

pdb_root *ssbuild(pdb_root *pdb, const char *ss_name)
{
  int serNum = 0;

  float dist;

  pdb_atom *curr_atom1, *curr_atom2;
  pdb_residue *curr_residue1, *curr_residue2;

  pdb_ssbond *ssbond;



  pdb->ssbonds = NULL;

  for (curr_atom1 = pdb->first_chain->first_residue->first_atom;
       curr_atom1->residue->next;
       curr_atom1 = curr_atom1->next) {

    curr_residue1 = curr_atom1->residue;

    if (NOT_CYS(curr_residue1, curr_atom1))
      continue;

    for (curr_atom2 = curr_atom1->residue->next->first_atom;
	 curr_atom2;
	 curr_atom2 = curr_atom2->next) {

      curr_residue2 = curr_atom2->residue;

      if (NOT_CYS(curr_residue2, curr_atom2))
	continue;

      dist = vecDist(curr_atom1->pos, curr_atom2->pos);

      if (dist < MAX_SSDIST) {
	strncpy(curr_residue1->resName, ss_name, PDB_RES_NAME_LEN-1);
	curr_residue1->resName[PDB_RES_NAME_LEN-1] = '\0';

	strncpy(curr_residue2->resName, ss_name, PDB_RES_NAME_LEN-1);
	curr_residue2->resName[PDB_RES_NAME_LEN-1] = '\0';

	serNum++;
	pdb->ssbonds = reallocate(pdb->ssbonds,
				  (serNum+1) * sizeof(*pdb->ssbonds));
	ssbond = allocate(sizeof(*ssbond));

	ssbond->serNum = serNum;

	ssbond->ss1.chainID = curr_residue1->chain->chainID;
	ssbond->ss1.seqNum = curr_residue1->resSeq;
	ssbond->ss1.icode = curr_residue1->iCode;

	ssbond->ss2.chainID = curr_residue2->chain->chainID;
	ssbond->ss2.seqNum = curr_residue2->resSeq;
	ssbond->ss2.icode = curr_residue2->iCode;

	ssbond->ss1.SymOP[0] = ssbond->ss2.SymOP[0] = '\0';

	ssbond->Length = (float) sqrt(dist);

	pdb->ssbonds[serNum-1] = ssbond;

	break;
      }
    }
  }

  if (pdb->ssbonds) {
    pdb->ssbonds[serNum] = NULL;
  }

  return pdb;    
}

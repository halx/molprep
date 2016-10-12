/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Driver program for the hbuild code.  Reads input file from first command line
 * argument or stdin.
 *
 *
 * $Id: molprep.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <errno.h>
#include <libgen.h>
#include <limits.h>

#include "common.h"
#include "pdb.h"
#include "top.h"
#include "ssbuild.h"
#include "protonate.h"
#include "hbuild.h"
#include "config.h"
#include "util/hashtab.h"
#include "util/util.h"


#define INPUT_LINE_LEN 258
#define PDB_TYPE_LEN 7
#define INPUT_DELIMITER  " =\t\n"
#define FILE_REQ(f,t) if (*f == '\0' )	\
    prerror(1, "%s input file required.\n", t)

#ifndef PATH_MAX
#define PATH_MAX 256		/* need a more portable constant */
#endif


#define X(a, b, c) c,
struct opt_flags options = {
#include "options.def"
};
#undef X

struct _opt_dict {
  const char *key;
  bool *val;
};



/*
 * scmp: comparison function for qsort(3) and bsearch(3), see C-FAQ 13.8 and 13.9
 *
 * in:  two general pointers
 * out: -1 if first argument is less than second, 0 if equal, 1 if larger
 *
 */

static int scmp(const void *p1, const void *p2)
{
  const struct _opt_dict *sp1 = (const struct _opt_dict *) p1;
  const struct _opt_dict *sp2 = (const struct _opt_dict *) p2;

  return strcmp(sp1->key, sp2->key);
}


int main(int argc, char **argv)
{
  char altloc_ind = 'A';

  int line_cnt = 0, model_no = INT_MIN, nssb = 0;

  char pdb_in_filename[PATH_MAX] = "\0";
  char pdb_out_filename[PATH_MAX] = "\0";
  char top_filename[PATH_MAX] = "\0";
  char ttb_filename[PATH_MAX] = "\0";
  char pdb_std_out_type[PDB_TYPE_LEN] = "\0";
  char ss_name[PDB_RES_NAME_LEN] = "CYS2";
  char buffer[INPUT_LINE_LEN];

  char *progname;
  char *key, *val, *bufp, *end;

  float pH = 7.0;

  FILE* input_stream;

  struct _opt_dict *od;

  pdb_root *pdb = NULL;
  topol_hash *top = NULL;

#define X(a, b, c) {a, b},
struct _opt_dict opt_dict[] = {
#include "options.def"
  {NULL, NULL}
};
#undef X

  const size_t opt_dict_size = sizeof(opt_dict) /
                               sizeof(struct _opt_dict) - 1;



  fprintf(stdout, "=== molprep %i.%i.%i ===\n\n", VERSION_MAJOR,
          VERSION_MINOR, VERSION_PATCH);

  progname = basename(argv[0]);

  if (argv[1]) {
    if (!(input_stream = fopen(argv[1], "r")) ) {
      perror(argv[1]);
      exit(EXIT_FAILURE);
    }
  } else {
    input_stream = stdin;
  }

  qsort(opt_dict, opt_dict_size, sizeof(struct _opt_dict), scmp);

  while (fgets(buffer, INPUT_LINE_LEN, input_stream) ) {
    line_cnt++;

    if ( !(bufp = normln(buffer) ) )
      continue;

    if ( !(key = strtok(bufp, INPUT_DELIMITER) ) ) {
      prerror(1, "%s: no key found in input (line %d).\n",
	      progname, line_cnt);
    }

    if ( !(val = strtok(NULL, INPUT_DELIMITER) ) ) {
      prerror(1, "%s: no value found in input (line %d).\n", 
	      progname, line_cnt);
    }

    if (STREQ(key, "inPDB") ) {
      strncpy(pdb_in_filename, val, PATH_MAX-1);
      pdb_in_filename[PATH_MAX-1] = '\0';
    } else if (STREQ(key, "outPDB") ) {
      strncpy(pdb_out_filename, val, PATH_MAX-1);
      pdb_out_filename[PATH_MAX-1] = '\0';
    } else if (STREQ(key, "top_file") ) {
      strncpy(top_filename, val, PATH_MAX-1);
      top_filename[PATH_MAX-1] = '\0';
    } else if (STREQ(key, "output_format") ) {
      strncpy(pdb_std_out_type, val, PDB_TYPE_LEN-1);
      pdb_std_out_type[PDB_TYPE_LEN-1] = '\0';
    } else if (STREQ(key, "altloc") ) {
      altloc_ind = *val;
    } else if (STREQ(key, "ss_name") ) {
	if (!pdb_format_residue(ss_name, val) )
	  prerror(1, "%s: ss_name cannot be longer than %d characters (line %d).\n",
		  progname, PDB_RES_NAME_LEN-1, line_cnt);
    } else if (STREQ(key, "model_no") ) {
      model_no = atoi(val);
    } else if (STREQ(key, "protonate_ttb") ) {
      strncpy(ttb_filename, val, PATH_MAX-1);
      ttb_filename[PATH_MAX-1] = '\0';
    } else if (STREQ(key, "protonate_pH") ) {
      errno = 0;
      pH = strtof(val, &end);

      if (end == val || errno == ERANGE) {
	prerror(1, "%s: cannot convert pH (line %d).\n", progname, line_cnt);
      }

      if (pH <= 0.0 || pH >= 14.0) {
	prwarn("extreme pH = %.2f.\n", pH);
      }
    } else {
      if((od = bsearch(&key, opt_dict, opt_dict_size,
		       sizeof(struct _opt_dict), scmp) ) ) {

	*val = tolower((unsigned char) *val);

	if (*val == 'y' || *val == 't')
	  *(od->val) = true;
	else
	  *(od->val) = false;
      } else {
	prerror(1, "%s, Unknown parameter in line %d: %s.\n",
		progname, line_cnt, key);
      }
    }
  }

  if (input_stream != stdin)
    fclose(input_stream);

  FILE_REQ(pdb_in_filename, "PDB input");
  FILE_REQ(pdb_out_filename, "PDB output");

  if (*top_filename == '\0') {
    strncpy(top_filename, TOP_DEFAULT_FILE, PATH_MAX - 1);
    top_filename[PATH_MAX-1] = '\0';
  }

  if (options.prot) {
    FILE_REQ(ttb_filename, "titratable translation table");
  }

  top = top_read(top, top_filename);
  pdb = pdb_read(pdb, pdb_in_filename, ss_name, model_no, &nssb);

  if (!options.rssb && nssb > 1)
    pdb = ssbuild(pdb, ss_name);

  if (options.prot)
    protonate(pdb, top->hash_table, ttb_filename, pH, altloc_ind);

  hbuild(pdb, top->hash_table, altloc_ind);

  pdb_write(pdb, pdb_out_filename, pdb_std_out_type, ss_name, altloc_ind);

  top_destroy(top);
  top = NULL;

  pdb_destroy(pdb);
  pdb = NULL;
}

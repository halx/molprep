/*
 * test interface to PROPKA 2.0
 *
 *
 * compile like:
 *
 * gfortran -Wall -Wextra -g3 -c propka.f 
 * gcc -g3 -Wall -Wextra -c test_f77_interface.c
 * gfortran -g3 -o test_f77_interface test_f77_interface.o propka.o
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libgen.h>
#include <unistd.h>

#include "propka.h"

#define STRNEQ(a,b,n) ( !strncmp((a), (b), (n)) )

#define INPUT_LINE_LEN 258
#define OUT_FILE "propka.out"
#define RETSTR_SIZE 65536


int main(int argc, char **argv)
{
  int natom;

  char *progname, *retstr;
  char buffer[INPUT_LINE_LEN];

  FILE* input_stream;



  progname = basename(argv[0]);

  if (argc < 2) {
    fprintf(stderr, "Usage: %s pdb_file\n", progname);
    exit(EXIT_FAILURE);
  }

  if (!(input_stream = fopen(argv[1], "r")) ) {
    perror(argv[1]);
    exit(EXIT_FAILURE);
  }

  natom = 0;

  while (fgets(buffer, INPUT_LINE_LEN, input_stream) ) {
    if (STRNEQ(buffer, "ATOM", 4) || STRNEQ(buffer, "HETATM", 6)) {
      natom++;
    }
  }

  retstr = malloc(RETSTR_SIZE);
  retstr[RETSTR_SIZE-1] = '\0';
  memset(retstr, ' ', RETSTR_SIZE);	/* propka may write NUL chars... */

  unlink("propka.out");

  runpka_(&natom, argv[1], OUT_FILE, retstr,
	  strlen(argv[1]), strlen(OUT_FILE), RETSTR_SIZE);
  printf("%s\n", retstr);

  return 0;
}

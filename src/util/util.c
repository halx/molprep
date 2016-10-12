/*
 * Copyright (C) 2011 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: util.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>

#include "common.h"
#include "vec.h"


#define COMMENT_CHAR '#'
#define QUOTE_CHAR '\\'



/*
 * simple fprintf(3) wrappers for error, warning and note message printing
 *
 */

void prerror(int status, const char* format, ... )
{
  va_list args;


  fprintf(stderr, "E> ");

  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);

  if (status) {
    exit(status);
  }
}

#define OUT stdout

void prwarn(const char* format, ... )
{
  va_list args;


  fprintf(OUT, "W> ");

  va_start(args, format);
  vfprintf(OUT, format, args);
  va_end(args);
}

void prnote(const char* format, ... )
{
  va_list args;


  fprintf(OUT, "N> ");

  va_start(args, format);
  vfprintf(OUT, format, args);
  va_end(args);
}

#ifndef NDEBUG
void dbgprintf(const char *funcname, const char *filename, unsigned int lineno,
	       const char *format, ...)
{
  va_list args;


  fprintf(stderr, "<%s() in %s, l.%i> ", funcname, filename, lineno);

  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
}
#endif


/*
 * simple malloc(3) and realloc(3) wrappers
 *
 */

void *allocate(size_t size)
{
  void *tmp;


  tmp = malloc(size);

  if (!tmp) {
    prerror(64, "Cannot allocate more memory (%d requested)\n", size);
  }

  return tmp;
}

void *reallocate(void *ptr, size_t size)
{
  void *tmp;


  tmp = realloc(ptr, size);

  if (!tmp) {
    prerror(65, "Cannot reallocate more memory (%d requested)\n", size);
  }

  return tmp;
}


/*
 * normln: normalise string by skipping leading white-space and
 *         removing comments and trailing white-space
 *
 * in:  null-terminated string
 * out: normalised string or NULL if null-string or NULL pointer passed
 *
 */

char *normln(char *string)
{
  int backslash, bkcnt;

  char *strp, *buffer, *bufp;


  if (!string || !*string) {
    return NULL;
  }

  /* skip leading white space */
  while (isspace((unsigned char) *string) && *string++)
    ;

  if (!*string) {
    return NULL;
  }

  /* delete comments */
  buffer = bufp = allocate(strlen(string) + 1);
  backslash = bkcnt = 0;

  for (strp = string; *strp; strp++) {
    if (!backslash) {
      if (*strp == QUOTE_CHAR) {
	backslash = 1;
	bkcnt++;
	continue;
      }

      if (*strp == COMMENT_CHAR) {
	*strp = '\0';
	break;
      }
    }

    backslash = 0;
    *bufp++ = *strp;
  }

  *bufp = '\0';

  if (bkcnt > 0) {
    strcpy(string, buffer);
    strp -= bkcnt;
  }

  free(buffer);

  /* delete trailing white space */
  while (isspace((unsigned char) *--strp) && strp >= string) {
    *strp = '\0';
  }

  if (!*string) {
    return NULL;
  } else {
    return string;
  }
}


/*
 * ishydrogen: check if atom is a hydrogen
 *
 * in:  element and atom name
 * out: true if hydrogen, false otherwise
 *
 */

int ishydrogen(const char *element, const char *name)
{
  while (isspace((unsigned char) *name) && *name++)
    ;

  return ((element[0] == ' ' && element[1] == 'H') ||
	  name[0] == 'H' ||
	  (isdigit((unsigned char) name[0]) && name[1] == 'H') );
}


/*
 * delnl: delete newline from an fgets read
 *
 * in: string
 *
 */

void delnl(char *string)
{
  char *s = string + strlen(string) - 1;

  if (*s == '\n') {
    *s = '\0';
  }
}


/*
 * itoa and reverse from K&R2
 */

/* reverse:  reverse string s in place */
void reverse(char s[])
{
  char c;
 
  for (size_t i = 0, j = strlen(s)-1; i<j; i++, j--) {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}
/* itoa and reverse from K&R2 - END */

/* itoa:  convert n to characters in s */
void itoa(int n, char s[])
{
  int i, sign;
 
  if ((sign = n) < 0)  /* record sign */
    n = -n;          /* make n positive */
  i = 0;
  do {       /* generate digits in reverse order */
    s[i++] = n % 10 + '0';   /* get next digit */
  } while ((n /= 10) > 0);     /* delete it */
  if (sign < 0)
    s[i++] = '-';
  s[i] = '\0';
  reverse(s);
}


/*
 * vector arithmetic helper functions
 *
 */

void vecCreate(fvec r, const float x, const float y, const float z)
{
  r[0] = x;
  r[1] = y;
  r[2] = z;
}

void vecCopy(fvec r, const fvec v)
{
  r[0] = v[0];
  r[1] = v[1];
  r[2] = v[2];
}

float vecLen(const fvec r)
{
  return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}

void vecScalarMult(fvec r, const float a, const fvec v)
{
  r[0] = a * v[0];
  r[1] = a * v[1];
  r[2] = a * v[2];
}

void vecScalarDiv(fvec r, const fvec v, const float a)
{
  r[0] = v[0] / a;
  r[1] = v[1] / a;
  r[2] = v[2] / a;
}

void vecAdd(fvec r, const fvec v, const fvec w)
{
  r[0] = v[0] + w[0];
  r[1] = v[1] + w[1];
  r[2] = v[2] + w[2];
}

void vecSub(fvec r, const fvec v, const fvec w)
{
  r[0] = v[0] - w[0];
  r[1] = v[1] - w[1];
  r[2] = v[2] - w[2];
}

float vecDist(const fvec c1, const fvec c2)
{
  float xd, yd, zd;

  xd = c1[0] - c2[0];
  yd = c1[1] - c2[1];
  zd = c1[2] - c2[2];

  return xd * xd + yd * yd + zd * zd;
}

void vecCrossProd(fvec r, const fvec v, const fvec w)
{
  r[0] = v[1] * w[2] - v[2] * w[1];
  r[1] = v[2] * w[0] - v[0] * w[2];
  r[2] = v[0] * w[1] - v[1] * w[0];
}

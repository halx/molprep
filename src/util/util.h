/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: util.h 161 2012-06-25 12:51:40Z hhl $
 *
 */


#ifndef _UTIL_H
#define _UTIL_H      1

#include "vec.h"

char *normln (char *string);
int ishydrogen(const char *element, const char *name);
void delnl(char *string);

void reverse(char s[]);
void itoa(int n, char s[]);

void vecCreate(fvec r, const float x, const float y, const float z);
void vecCopy(fvec r, const fvec v);
float vecLen(const fvec r);
void vecScalarMult(fvec r, const float a, const fvec v);
void vecScalarDiv(fvec r, const fvec v, const float a);
void vecAdd(fvec r, const fvec v, const fvec w);
void vecSub(fvec r, const fvec v, const fvec w);
float vecDist(const fvec c1, const fvec c2);
void vecCrossProd(fvec r, const fvec v, const fvec w);

void *allocate(size_t size);
void *reallocate(void *ptr, size_t size);

void prerror(int status, const char* format, ... );
void prwarn(const char* format, ... );
void prnote(const char* format, ... );

#endif

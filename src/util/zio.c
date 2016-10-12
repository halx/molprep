/*
 * Copyright (C) 2011 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Zlib I/O support via wrapper functions.
 *
 *
 * $Id: zio.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include "config.h"

#include <stdio.h>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif



void *fzopen(const char *path, const char *mode)
{
  if (!path) return NULL;

#ifdef HAVE_ZLIB
  return gzopen(path, mode);
#else
  return fopen(path, mode);
#endif
}

char *fzgets(void *file, char *s, int size)
{
  if (!file || !s || size < 3) return NULL;

#ifdef HAVE_ZLIB
  return gzgets(file, s, size);
#else
  return fgets(s, size, file);
#endif
}

int fzclose(void *file)
{
#ifdef HAVE_ZLIB
  return gzclose(file);
#else
  return fclose(file);
#endif
}

/*
 * Copyright (C) 2011 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Header file for Zlib I/O support.
 *
 *
 * $Id: zio.h 161 2012-06-25 12:51:40Z hhl $
 *
 */



void *fzopen(const char *path, const char *mode);
char *fzgets(void *file, char *s, int size);	 /* parameter order as gzgets */
int fzclose(void *file);

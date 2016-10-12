/*
 * Copyright (C) 2011-2012 Hannes Loeffler, STFC Daresbury, UK
 *
 * $Id: common.h 153 2012-06-18 10:18:32Z hhl $
 *
 */



#ifndef _COMMON_H
#define _COMMON_H      1

#include <stdio.h>
#include <string.h>
#include <stdbool.h>


#ifndef NDEBUG
void dbgprintf(const char *funcname, const char *filename, unsigned int lineno,
	       const char *format, ...);
#define DEBUG_PRINT(...) do {						\
    dbgprintf(__func__, __FILE__, __LINE__, __VA_ARGS__); } while (0)
#else
#define DEBUG_PRINT(...) do { } while (0)
#endif

#define STREQ(a,b) ( !strcmp((a), (b)) )
#define STRNEQ(a,b,n) ( !strncmp((a), (b), (n)) )


/* global structure to hold the option flags */
extern struct opt_flags {
  bool remh, nomodel, nocryst, noter, noend, prot, rssb, wrss, keepssn,
    keepser, nterm, cterm, dna5term, dna3term, rna5term, rna3term, warnocc;
} options;

#endif

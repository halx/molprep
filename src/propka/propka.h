/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Interface to PROPKA 2.00 (2008-11-12).
 *
 *
 * $Id: propka.h 150 2012-06-14 09:56:02Z hhl $
 *
 */



#ifndef _PROPKA_H
#define _PROPKA_H      1

int runpka_(unsigned int *maxatm, unsigned int *maxres, unsigned int *maxar,
	    char *pdbfile, char *outfile, char *retstr,
	    unsigned int pdb_len, unsigned int out_len, unsigned int r_len);

#endif	/* !_PROPKA_H */

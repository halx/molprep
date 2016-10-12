/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal dynamic array implementation.
 *
 * NOTE: it is the caller's resonsibility to allocate memory for the data.
 *
 *
 * $Id: darray.h 161 2012-06-25 12:51:40Z hhl $
 *
 */



#ifndef _DARRAY_H
#define _DARRAY_H      1


#define FOREACH_DARRAY(i, darray)		\
  for (unsigned int i = 0; i < darray_num_el(darray); i++)

typedef struct _DArray DArray;

DArray *darray_init(DArray *darray, size_t initrate);
void darray_add(DArray *darray, void *data);
void *darray_get(DArray *darray, unsigned int idx);
size_t darray_num_el(DArray *darray);
void *darray_shrink(DArray *darray);
void darray_destroy(DArray *darray);

#endif

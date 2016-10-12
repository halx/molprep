/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal dynamic array implementation.
 *
 * NOTE: it is the caller's resonsibility to allocate memory for the data.
 *
 *
 * $Id: darray.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <stdlib.h>

#include "common.h"
#include "darray.h"
#include "util.h"



struct _DArray {
  void **data;			// pointer to data
  size_t nel;			// number of elements
  size_t max;			// maximum number of elements
  size_t exrate;		// expansion rate
};


DArray *darray_init(DArray *darray, size_t initmax)
{
  if (darray)
    return darray;

  darray = allocate(initmax * sizeof *darray);
  darray->data = allocate(initmax * sizeof(void *) );

  darray->nel = 0;
  darray->max = initmax;
  darray->exrate = initmax;	// linear expansion rate

  return darray;
}


static void darray_expand(DArray *darray, size_t newsize)
{
  size_t old_max = darray->max; 


  darray->data = reallocate(darray->data, newsize * sizeof(void *) );
  darray->max = newsize;

  memset(darray->data + old_max, 0, darray->exrate + 1);
}


// data is not allocated, only pointed to
void darray_add(DArray *darray, void *data)
{

  if (!darray || !data)
    return;

  if (darray->nel >= darray->max) {
    darray_expand(darray, darray->max + darray->exrate);
  }

  darray->data[darray->nel] = data;
  darray->nel++;
}


void *darray_get(DArray *darray, unsigned int idx)
{
  if (!darray || !darray->data || idx > darray->max)
    return NULL;
  else
    return darray->data[idx];
}


size_t darray_num_el(DArray *darray)
{
  if (darray)
    return darray->nel;
  else
    return 0;
}


void darray_destroy(DArray *darray)
{
  if (darray) {
    if (darray->data)
      free(darray->data);

    free(darray);
  }
}

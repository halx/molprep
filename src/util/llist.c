/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal linked list implementation.
 *
 * NOTE: it is the caller's resonsibility to allocate memory for the data.
 *
 *
 * $Id$
 *
 */



#include <stdlib.h>

#include "common.h"
#include "llist.h"
#include "util.h"


struct _Listnode {
  void *restrict data;
  struct _Listnode *restrict prev;
  struct _Listnode *restrict next;
};

struct _List {
  unsigned int nel;
  Listnode *restrict first;
  Listnode *restrict last;
};


List *list_init(List *restrict list)
{
  if (list)
    return list;

  list = allocate(sizeof *list);

  list->nel = 0;
  list->first = list->last = NULL;

  return list;
}


static List *add_node(List *restrict list, void *restrict data)
{
  Listnode *new;


  new = allocate(sizeof new);
  new->data = data;
  list->nel++;

  if (!list->last) {
    list->first = new;
    new->prev = NULL;
  } else {
    list->last->next = new;
    new->prev = list->last; 
  }

  new->next = NULL;
  list->last = new;

  return list;
}


void list_add(List *restrict list, void *restrict data)
{
  if (!list || !data)
    return;

  list = add_node(list, data);
}


// keeps only unique elements like a set
// first/last must be NULL on first invocation!
void list_add_uniq(List *restrict list, void *restrict data, size_t len)
{
  Listnode *new;


  if (!list || !data || len < 1)
    return;

  for (new = list->first; new; new = new->next) {
    // probably not efficient for non-byte data!
    if (STRNEQ(new->data, data, len)) {
      return;
    }
  }

  list = add_node(list, data);
}


Listnode *list_first(List *restrict list)
{
  if (!list)
    return NULL;
  else
    return list->first;
}

Listnode *list_next(Listnode *restrict node)
{
  if (!node)
    return NULL;
  else
    return node->next;
}


void *list_get(Listnode *restrict node)
{
  if (!node || !node->data)
    return NULL;
  else
    return node->data;
}


void *list_pop_first(List *restrict list)
{
  void *data;
  Listnode *tmp;


  if (!list || !list->first)
    return NULL;

  data = list->first->data;
  tmp = list->first->next;
  free(list->first);
  list->first = tmp;

  if (NULL == list->first)
    list->last = list->first;

  return data;
}


void *list_pop_last(List *restrict list)
{
  void *data;
  Listnode *tmp;


  if (!list || !list->last)
    return NULL;

  data = list->last->data;
  tmp = list->last->prev;
  free(list->last);
  list->last = tmp;

  if (NULL == list->last)
    list->first = list->last;

  return data;
}


unsigned int list_num_el(List *restrict list)
{
  if (!list || !list->first || !list->last)
    return 0;
  else
    return list->nel;
}


unsigned int list_is_empty(List *restrict list)
{
  if (!list || !list->first || !list->last) {
    return 1;
  } else {
    return 0;
  }
}


void list_destroy(List *restrict list)
{
  Listnode *node, *tmp;


  if (!list)
    return;

  node = list->first;

  while (node) {
    tmp = node;
    node = node->next;
    free(tmp);
  }

  free(list);
}

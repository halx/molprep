/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal linked list implementation.
 *
 *
 * $Id$
 *
 */



#ifndef _LLIST_H
#define _LLIST_H      1


#define FOREACH_LIST(node, list)					\
  for (node = list_first(list); node; node = list_next(node))

typedef struct _Listnode Listnode;
typedef struct _List List;

List *list_init(List *restrict list);
void list_add(List *restrict list, void *restrict data);
void list_add_uniq(List *restrict q, void *restrict data, size_t len);
Listnode *list_first(List *restrict list);
Listnode *list_next(Listnode *restrict node);
void *list_get(Listnode *restrict node);
void *list_pop_first(List *restrict list);
void *list_pop_last(List *restrict list);
unsigned int list_num_el(List *restrict list);
unsigned int list_is_empty(List *restrict list);
void list_destroy(List *restrict list);

#endif

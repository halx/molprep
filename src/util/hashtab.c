/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Hash table algorithms.  Only fixed table size possible at the moment.
 *
 *
 * $Id: hashtab.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <string.h>
#include <stdlib.h>

#include "common.h"
#include "hashtab.h"
#include "hashfuncs.h"
#include "util.h"


struct _Hashnode {
  char *restrict key;
  void *restrict data;
  struct _Hashnode *restrict next;
};

struct _Hashtable {
  Hashnode **entries;
  hashfp hash_func;
  unsigned int table_size;
};



Hashtable *hash_init(const hashfp hash_func, const unsigned int table_size)
{
  Hashtable *table;


  table = allocate(sizeof(*table));
  table->entries = allocate(table_size * sizeof(Hashnode));

  for (unsigned int i = 0; i < table_size; i++) {
    table->entries[i] = NULL;
  }

  table->hash_func = hash_func;
  table->table_size = table_size;

  return table;
}

/* FIXME: can hash be negative? */
void hash_insert(const Hashtable *restrict table, char *restrict key,
		 const size_t key_len, void *restrict data)
{
  unsigned int idx;
  Hashnode *new, *old = NULL;


  idx = table->hash_func(key) & (table->table_size - 1);

  /* silently ignore keys already existing */
  for (Hashnode *tp = table->entries[idx]; tp; old = tp, tp = tp->next) {
    if(STRNEQ(tp->key, key, key_len - 1) ) {
      return;
    }
  }

  new = allocate(sizeof(*new) );
  new->key = key;
  new->data = data;
  new->next = NULL;

  if (table->entries[idx]) {
    old->next = new;		/* append */
  } else {
    table->entries[idx] = new;
  }
}

Hashnode *hash_search(const Hashtable *table, const char *key,
		       const size_t key_len)
{
  unsigned int idx;


  idx = table->hash_func(key) & (table->table_size - 1);

  for (Hashnode *tp = table->entries[idx]; tp; tp = tp->next) {
    if(STRNEQ(key, tp->key, key_len - 1) ) {
      return tp;
    }
  }
 
  return NULL;
}

#ifndef NDEBUG
#include <stdio.h>

void hash_print(const Hashtable *table)
{

  for (unsigned int i = 0; i < table->table_size; i++) {
    printf("%i:", i);

    for (Hashnode *tp = table->entries[i]; tp; tp = tp->next) {
      printf(" %s", tp->key);
    }

    printf("\n");
  }
}
#endif

/* getter to support encapsulation */
void *hash_node_get_data(const Hashnode *node)
{
  return node->data;
}

void hash_destroy(Hashtable *table)
{
  Hashnode *tp, *old;


  for (unsigned int i = 0; i < table->table_size; i++) {
    tp = table->entries[i];

    while(tp) {
      old = tp;
      tp = tp->next;
      free(old);
    }
  }

  free(table->entries);
  free(table);
}

/* determin high bit for correct table size */
int hibit(unsigned int n) {

  n |= (n >>  1);
  n |= (n >>  2);
  n |= (n >>  4);
  n |= (n >>  8);
  n |= (n >> 16);

  return n - (n >> 1);
}

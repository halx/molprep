/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Header file for hash function algorithms.
 *
 *
 * $Id: hashtab.h 161 2012-06-25 12:51:40Z hhl $
 *
 */



#ifndef _HASH_H
#define _HASH_H      1

#include <stdint.h>


typedef uint32_t (*hashfp)(const char *str);
typedef struct _Hashnode Hashnode;
typedef struct _Hashtable Hashtable;

Hashtable *hash_init(const hashfp hash_func, const unsigned int table_size);
void hash_insert(const Hashtable *restrict table, char *restrict key,
		 const size_t key_len, void *restrict data);
Hashnode *hash_search(const Hashtable *table, const char *key,
		       const size_t key_len);

#ifndef NDEBUG
void hash_print(const Hashtable *table);
#endif

void *hash_node_get_data(const Hashnode *node);
void hash_destroy(Hashtable *table);
int hibit(unsigned int n);

#endif

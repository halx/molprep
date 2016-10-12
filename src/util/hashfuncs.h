/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * Header file for library of simple hash functions.
 *
 *
 * $Id: hashfuncs.h 161 2012-06-25 12:51:40Z hhl $
 */


#ifndef _HASHFUNCS_H
#define _HASHFUNCS_H      1

#include <stdint.h>

#ifndef UINT32_MAX
#define uint32_t uint_least32_t
#endif

#ifndef UINT8_MAX
#define uint8_t uint_least8_t
#endif


uint32_t kandr2_hash(const char *str);
uint32_t pl_hash(const char *str);
uint32_t djb2_hash(const char *str);
uint32_t rs_hash(const char *str);
uint32_t sdbm_hash(const char *str);
uint32_t ap_hash(const char *str);
uint32_t fnv1a_hash(const char *str);
uint32_t oat_hash(const char *str);
uint32_t sbox_hash(const char *str);
uint32_t ly_hash(const char *str);
uint32_t am_hash(const char *str);
uint32_t rot13_hash(const char *str);
uint32_t crc32_hash(const char *str);
uint32_t dek_hash(const char *str);
uint32_t lookup3_hash(const char *str);
uint32_t murmur2_hash(const char *str);
uint32_t ph_hash(const char *str);
uint32_t dh_hash(const char *str);
uint32_t jesteress_hash(const char *str);
uint32_t meiyan_hash(const char *str);
uint32_t mantis_hash(const char *str);
uint32_t whiz_hash(const char *str);
uint32_t alfalfa_hash(const char *str);

#endif

/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A library of simple hash functions.  These functions work (reasonably) well
 * with the 4 character PDB keys ([ 0-9A-Z]) but functions are dependant on the
 * input, i.e. a different function may be somewhat more efficient if the set
 * of keys changes.  Other functions may be required for longer keys, see e.g.
 * http://www.strchr.com/hash_functions .
 *
 * All functions have the same parameters.  The string (first parameter) is
 * required to be \0 terminated.  The table size (second parameter) must be 2^n.
 * The returned hash value must be suitably mod'ed to fit the required table
 * size.
 *
 *
 * $Id: hashfuncs.c 161 2012-06-25 12:51:40Z hhl $
 *
 */



#include <stdint.h>

#include "hashfuncs.h"


/* K&R 2; the book uses 31 as constant but 13331 seems to work better here */
uint32_t kandr2_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while ( (c = *str++) )
    hash = 13331 * hash + c;

  return hash;
}

/* Paul Larson; like K&R 2 but uses 101 as constant */
uint32_t pl_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while ( (c = *str++) )
    hash = 101 * hash + c;

  return hash;
}

/* Ozan Yigit
   from sdbm (a public-domain reimplementation of ndbm) database library.*/
uint32_t sdbm_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while ( (c = *str++) )
    hash = (hash << 6) + (hash << 16) - hash + c;  /* 65599 * hash + c */

  return hash;
}

/* Daniel J. Bernstein 2; + instead of ^ in previous version
   (actually from Gosling Emacs mid-1980s???)
   http://cr.yp.to/djbdns.html */
uint32_t djb2_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 5381U;


  while ( (c = *str++) )
    hash = ((hash << 5) + hash) ^ c;  /* hash = 33 * hash ^ c */

  return hash;
}

/* Leonid Yuriev */
uint32_t ly_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;
 

  while( (c = *str++) )
    hash = (hash * 1664525U) + c + 1013904223U;

  return hash;
}

/* from Robert Sedgwick's "Algorithms in C" book */
uint32_t rs_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;
  /* uint32_t a = 31415U, b = 27183U;  // "universial" hash */
  uint32_t a = 63689U, b = 378551U;


  while( (c = *str++) ) {
    hash = hash * a + c;
    a *= b;
  }

  return hash;
}

/* Serge Vakulenko
   http://vak.ru/doku.php/proj/hash/efficiency-en */
uint32_t rot13_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while( (c = *str++) ) {
    hash += c;
    hash -= (hash << 13) | (hash >> 19);
  }

  return hash;
}

/* Arash Partow
   Free use  permitted under the guidelines and in accordance with the most
   current version of the "Common Public License",
   http://www.opensource.org/licenses/cpl1.0.php */
uint32_t ap_hash(const char *str)
{
  uint8_t c;
  uint32_t i = 0, hash = 0;


  while( (c = *str++) ) {
    i++;

    hash ^= ((i & 1) == 0) ?
      ( (hash << 7) ^ c ^ (hash >> 3)) :
      (~((hash << 11) ^ c ^ (hash >> 5)));
  }

  return hash;
}

/* Fowler, Noll, Vo's FNV-1a; in FNV-1 order is reversed
   in the public domain
   http://isthe.com/chongo/tech/comp/fnv/ */
uint32_t fnv1a_hash(const char *str)
{
  uint8_t c;
  uint32_t prime = 16777619U;
  uint32_t hash = 2166136261U;


  while( (c = *str++) )
    hash = (hash ^ c) * prime;

  return hash;
}

#include <string.h>

/* Donald E. Knuth */
uint32_t dek_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = strlen(str);

  while( (c = *str++) )
    hash = ((hash << 5) ^ (hash >> 27)) ^ c;

  return hash;
}

/* Bob Jenkin's One-at-a-time
   in the public domain
   http://burtleburtle.net/bob/hash/doobs.html */
uint32_t oat_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while ( (c = *str++) ) {
    hash += c;
    hash += ( hash << 10 );
    hash ^= ( hash >> 6 );
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  return hash;
}

/* Murmur2
   in the public domain for non-business purposes, otherwise MIT license
   https://code.google.com/p/smhasher/ */
uint32_t murmur2_hash(const char *str)
{
  /* 'm' and 'r' are mixing constants generated offline.
     They're not really 'magic', they just happen to work well. */
  const uint32_t m = 0x5bd1e995;
  const int r = 24;

  /* Initialize the hash to a 'random' value */
  uint32_t k, seed = 0x3FB0BB5F;
  size_t len = strlen(str);
  uint32_t h = seed ^ len;

  /* Mix 4 bytes at a time into the hash */
  const uint8_t *data = (const uint8_t *)str;


  while(len >= 4) {
    k = *(const uint32_t *)data;

    k *= m; 
    k ^= k >> r; 
    k *= m; 
		
    h *= m; 
    h ^= k;

    data += 4;
    len -= 4;
  }
	
  /* Handle the last few bytes of the input array */
  switch(len) {
  case 3: h ^= data[2] << 16;
  case 2: h ^= data[1] << 8;
  case 1: h ^= data[0];
    h *= m;
  };

  /* Do a few final mixes of the hash to ensure the last few
     bytes are well-incorporated. */
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;

  return h;
}

#define WORD uint16_t

/* Paul Hsieh
   LGPL 2.1 license
   http://www.azillionmonkeys.com/qed/hash.html */
uint32_t ph_hash(const char *str)
{
  size_t len = strlen(str);
  uint32_t tmp, rem, hash = len;


  if(len == 0) return 0;

  rem = len & 3;
  len >>= 2;

  /* Main loop */
  for (;len > 0; len--) {
    hash += *(const WORD *)str;
    tmp = (*(const WORD *) (str + 2) << 11) ^ hash;
    hash = (hash << 16) ^ tmp;
    str += 2 * sizeof (WORD);
    hash += hash >> 11;
  }

  /* Handle end cases */
  switch(rem) {
  case 3:
    hash += *(const WORD *)str;
    hash ^= hash << 16;
    hash ^= str[sizeof (WORD)] << 18;
    hash += hash >> 11;
    break;

  case 2:
    hash += *(const WORD *)str;
    hash ^= hash << 11;
    hash += hash >> 17;
    break;

  case 1:
    hash += *str;
    hash ^= hash << 10;
    hash += hash >> 1;
  }

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}

/* Bret Mulvey's SBox
   http://home.comcast.net/~bretm/hash/10.html */
static const uint32_t sbox_table[] = {
  0xf53e1837, 0x5f14c86b, 0x9ee3964c, 0xfa796d53, 0x32223fc3, 0x4d82bc98,
  0xa0c7fa62, 0x63e2c982, 0x24994a5b, 0x1ece7bee, 0x292b38ef, 0xd5cd4e56,
  0x514f4303, 0x7be12b83, 0x7192f195, 0x82dc7300, 0x084380b4, 0x480b55d3,
  0x5f430471, 0x13f75991, 0x3f9cf22c, 0x2fe0907a, 0xfd8e1e69, 0x7b1d5de8,
  0xd575a85c, 0xad01c50a, 0x7ee00737, 0x3ce981e8, 0x0e447efa, 0x23089dd6,
  0xb59f149f, 0x13600ec7, 0xe802c8e6, 0x670921e4, 0x7207eff0, 0xe74761b0,
  0x69035234, 0xbfa40f19, 0xf63651a0, 0x29e64c26, 0x1f98cca7, 0xd957007e,
  0xe71ddc75, 0x3e729595, 0x7580b7cc, 0xd7faf60b, 0x92484323, 0xa44113eb,
  0xe4cbde08, 0x346827c9, 0x3cf32afa, 0x0b29bcf1, 0x6e29f7df, 0xb01e71cb,
  0x3bfbc0d1, 0x62edc5b8, 0xb7de789a, 0xa4748ec9, 0xe17a4c4f, 0x67e5bd03,
  0xf3b33d1a, 0x97d8d3e9, 0x09121bc0, 0x347b2d2c, 0x79a1913c, 0x504172de,
  0x7f1f8483, 0x13ac3cf6, 0x7a2094db, 0xc778fa12, 0xadf7469f, 0x21786b7b,
  0x71a445d0, 0xa8896c1b, 0x656f62fb, 0x83a059b3, 0x972dfe6e, 0x4122000c,
  0x97d9da19, 0x17d5947b, 0xb1affd0c, 0x6ef83b97, 0xaf7f780b, 0x4613138a,
  0x7c3e73a6, 0xcf15e03d, 0x41576322, 0x672df292, 0xb658588d, 0x33ebefa9,
  0x938cbf06, 0x06b67381, 0x07f192c6, 0x2bda5855, 0x348ee0e8, 0x19dbb6e3,
  0x3222184b, 0xb69d5dba, 0x7e760b88, 0xaf4d8154, 0x007a51ad, 0x35112500,
  0xc9cd2d7d, 0x4f4fb761, 0x694772e3, 0x694c8351, 0x4a7e3af5, 0x67d65ce1,
  0x9287de92, 0x2518db3c, 0x8cb4ec06, 0xd154d38f, 0xe19a26bb, 0x295ee439,
  0xc50a1104, 0x2153c6a7, 0x82366656, 0x0713bc2f, 0x6462215a, 0x21d9bfce,
  0xba8eace6, 0xae2df4c1, 0x2a8d5e80, 0x3f7e52d1, 0x29359399, 0xfea1d19c,
  0x18879313, 0x455afa81, 0xfadfe838, 0x62609838, 0xd1028839, 0x0736e92f,
  0x3bca22a3, 0x1485b08a, 0x2da7900b, 0x852c156d, 0xe8f24803, 0x00078472,
  0x13f0d332, 0x2acfd0cf, 0x5f747f5c, 0x87bb1e2f, 0xa7efcb63, 0x23f432f0,
  0xe6ce7c5c, 0x1f954ef6, 0xb609c91b, 0x3b4571bf, 0xeed17dc0, 0xe556cda0,
  0xa7846a8d, 0xff105f94, 0x52b7ccde, 0x0e33e801, 0x664455ea, 0xf2c70414,
  0x73e7b486, 0x8f830661, 0x8b59e826, 0xbb8aedca, 0xf3d70ab9, 0xd739f2b9,
  0x4a04c34a, 0x88d0f089, 0xe02191a2, 0xd89d9c78, 0x192c2749, 0xfc43a78f,
  0x0aac88cb, 0x9438d42d, 0x9e280f7a, 0x36063802, 0x38e8d018, 0x1c42a9cb,
  0x92aaff6c, 0xa24820c5, 0x007f077f, 0xce5bc543, 0x69668d58, 0x10d6ff74,
  0xbe00f621, 0x21300bbe, 0x2e9e8f46, 0x5acea629, 0xfa1f86c7, 0x52f206b8,
  0x3edf1a75, 0x6da8d843, 0xcf719928, 0x73e3891f, 0xb4b95dd6, 0xb2a42d27,
  0xeda20bbf, 0x1a58dbdf, 0xa449ad03, 0x6ddef22b, 0x900531e6, 0x3d3bff35,
  0x5b24aba2, 0x472b3e4c, 0x387f2d75, 0x4d8dba36, 0x71cb5641, 0xe3473f3f,
  0xf6cd4b7f, 0xbf7d1428, 0x344b64d0, 0xc5cdfcb6, 0xfe2e0182, 0x2c37a673,
  0xde4eb7a3, 0x63fdc933, 0x01dc4063, 0x611f3571, 0xd167bfaf, 0x4496596f,
  0x3dee0689, 0xd8704910, 0x7052a114, 0x068c9ec5, 0x75d0e766, 0x4d54cc20,
  0xb44ecde2, 0x4abc653e, 0x2c550a21, 0x1a52c0db, 0xcfed03d0, 0x119bafe2,
  0x876a6133, 0xbc232088, 0x435ba1b2, 0xae99bbfa, 0xbb4f08e4, 0xa62b5f49,
  0x1da4b695, 0x336b84de, 0xdc813d31, 0x00c134fb, 0x397a98e6, 0x151f0e64,
  0xd9eb3e69, 0xd3c7df60, 0xd2f2c336, 0x2ddd067b, 0xbd122835, 0xb0b3bd3a,
  0xb0d54e46, 0x8641f1e4, 0xa0b38f96, 0x51d39199, 0x37a6ad75, 0xdf84ee41,
  0x3c034cba, 0xacda62fc, 0x11923b8b, 0x45ef170a
};

uint32_t sbox_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while ( (c = *str++) ) {
    hash ^= sbox_table[c];
    hash *= 3;
  }

  return hash;
}

/* CRC32
   http://www.faqs.org/rfcs/rfc1952.html */
static const uint32_t crc32_table[] = {
  0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f,
  0xe963a535, 0x9e6495a3, 0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
  0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91, 0x1db71064, 0x6ab020f2,
  0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
  0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9,
  0xfa0f3d63, 0x8d080df5, 0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
  0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b, 0x35b5a8fa, 0x42b2986c,
  0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
  0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423,
  0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
  0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d, 0x76dc4190, 0x01db7106,
  0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
  0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d,
  0x91646c97, 0xe6635c01, 0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
  0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950,
  0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
  0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7,
  0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
  0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa,
  0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
  0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81,
  0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
  0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683, 0xe3630b12, 0x94643b84,
  0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
  0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb,
  0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
  0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5, 0xd6d6a3e8, 0xa1d1937e,
  0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
  0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55,
  0x316e8eef, 0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
  0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28,
  0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
  0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f,
  0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
  0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242,
  0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
  0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69,
  0x616bffd3, 0x166ccf45, 0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
  0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc,
  0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
  0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693,
  0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
  0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
};

uint32_t crc32_hash(const char *str)
{
  uint8_t c;
  uint32_t crc = 0xffffffff;

 
  while ( (c = *str++) )
    crc = crc32_table[(crc ^ c) & 0xff] ^ (crc >> 8);

  return (crc ^ 0xffffffff);
}

/* Alexander Myasnikov, http://amsoftware.narod.ru/algo.html */
static const uint8_t am_table[] = {
  0xa3, 0xd7, 0x09, 0x83, 0xf8, 0x48, 0xf6, 0xf4, 0xb3, 0x21, 0x15, 0x78, 0x99,
  0xb1, 0xaf, 0xf9, 0xe7, 0x2d, 0x4d, 0x8a, 0xce, 0x4c, 0xca, 0x2e, 0x52, 0x95,
  0xd9, 0x1e, 0x4e, 0x38, 0x44, 0x28, 0x0a, 0xdf, 0x02, 0xa0, 0x17, 0xf1, 0x60,
  0x68, 0x12, 0xb7, 0x7a, 0xc3, 0xe9, 0xfa, 0x3d, 0x53, 0x96, 0x84, 0x6b, 0xba,
  0xf2, 0x63, 0x9a, 0x19, 0x7c, 0xae, 0xe5, 0xf5, 0xf7, 0x16, 0x6a, 0xa2, 0x39,
  0xb6, 0x7b, 0x0f, 0xc1, 0x93, 0x81, 0x1b, 0xee, 0xb4, 0x1a, 0xea, 0xd0, 0x91,
  0x2f, 0xb8, 0x55, 0xb9, 0xda, 0x85, 0x3f, 0x41, 0xbf, 0xe0, 0x5a, 0x58, 0x80,
  0x5f, 0x66, 0x0b, 0xd8, 0x90, 0x35, 0xd5, 0xc0, 0xa7, 0x33, 0x06, 0x65, 0x69,
  0x45, 0x00, 0x94, 0x56, 0x6d, 0x98, 0x9b, 0x76, 0x97, 0xfc, 0xb2, 0xc2, 0xb0,
  0xfe, 0xdb, 0x20, 0xe1, 0xeb, 0xd6, 0xe4, 0xdd, 0x47, 0x4a, 0x1d, 0x42, 0xed,
  0x9e, 0x6e, 0x49, 0x3c, 0xcd, 0x43, 0x27, 0xd2, 0x07, 0xd4, 0xde, 0xc7, 0x67,
  0x18, 0x89, 0xcb, 0x30, 0x1f, 0x8d, 0xc6, 0x8f, 0xaa, 0xc8, 0x74, 0xdc, 0xc9,
  0x5d, 0x5c, 0x31, 0xa4, 0x70, 0x88, 0x61, 0x2c, 0x9f, 0x0d, 0x2b, 0x87, 0x50,
  0x82, 0x54, 0x64, 0x26, 0x7d, 0x03, 0x40, 0x34, 0x4b, 0x1c, 0x73, 0xd1, 0xc4,
  0xfd, 0x3b, 0xcc, 0xfb, 0x7f, 0xab, 0xe6, 0x3e, 0x5b, 0xa5, 0xad, 0x04, 0x23,
  0x9c, 0x14, 0x51, 0x22, 0xf0, 0x29, 0x79, 0x71, 0x7e, 0xff, 0x8c, 0x0e, 0xe2,
  0x0c, 0xef, 0xbc, 0x72, 0x75, 0x6f, 0x37, 0xa1, 0xec, 0xd3, 0x8e, 0x62, 0x8b,
  0x86, 0x10, 0xe8, 0x08, 0x77, 0x11, 0xbe, 0x92, 0x4f, 0x24, 0xc5, 0x32, 0x36,
  0x9d, 0xcf, 0xf3, 0xa6, 0xbb, 0xac, 0x5e, 0x6c, 0xa9, 0x13, 0x57, 0x25, 0xb5,
  0xe3, 0xbd, 0xa8, 0x3a, 0x01, 0x05, 0x59, 0x2a, 0x46
};

uint32_t am_hash(const char *str)
{
  const uint32_t prime = 1717;
  const size_t len = strlen(str);
  uint32_t hash = len;

  for (uint32_t i = 0; i < len; i++) {
    hash ^= am_table[(str[i] + i) & 255];
    hash = hash * prime;
  }

  return hash;
}

/* David R. Hanson
   http://www.cs.princeton.edu/software/cii/doc/atom.html */
static const uint32_t dh_table[] = {
  2078917053, 143302914, 1027100827, 1953210302, 755253631, 2002600785,
  1405390230, 45248011, 1099951567, 433832350, 2018585307, 438263339,
  813528929, 1703199216, 618906479, 573714703, 766270699, 275680090,
  1510320440, 1583583926, 1723401032, 1965443329, 1098183682, 1636505764,
  980071615, 1011597961, 643279273, 1315461275, 157584038, 1069844923,
  471560540, 89017443, 1213147837, 1498661368, 2042227746, 1968401469,
  1353778505, 1300134328, 2013649480, 306246424, 1733966678, 1884751139,
  744509763, 400011959, 1440466707, 1363416242, 973726663, 59253759,
  1639096332, 336563455, 1642837685, 1215013716, 154523136, 593537720,
  704035832, 1134594751, 1605135681, 1347315106, 302572379, 1762719719,
  269676381, 774132919, 1851737163, 1482824219, 125310639, 1746481261,
  1303742040, 1479089144, 899131941, 1169907872, 1785335569, 485614972,
  907175364, 382361684, 885626931, 200158423, 1745777927, 1859353594,
  259412182, 1237390611, 48433401, 1902249868, 304920680, 202956538,
  348303940, 1008956512, 1337551289, 1953439621, 208787970, 1640123668,
  1568675693, 478464352, 266772940, 1272929208, 1961288571, 392083579,
  871926821, 1117546963, 1871172724, 1771058762, 139971187, 1509024645,
  109190086, 1047146551, 1891386329, 994817018, 1247304975, 1489680608,
  706686964, 1506717157, 579587572, 755120366, 1261483377, 884508252,
  958076904, 1609787317, 1893464764, 148144545, 1415743291, 2102252735,
  1788268214, 836935336, 433233439, 2055041154, 2109864544, 247038362,
  299641085, 834307717, 1364585325, 23330161, 457882831, 1504556512,
  1532354806, 567072918, 404219416, 1276257488, 1561889936, 1651524391,
  618454448, 121093252, 1010757900, 1198042020, 876213618, 124757630,
  2082550272, 1834290522, 1734544947, 1828531389, 1982435068, 1002804590,
  1783300476, 1623219634, 1839739926, 69050267, 1530777140, 1802120822,
  316088629, 1830418225, 488944891, 1680673954, 1853748387, 946827723,
  1037746818, 1238619545, 1513900641, 1441966234, 367393385, 928306929,
  946006977, 985847834, 1049400181, 1956764878, 36406206, 1925613800,
  2081522508, 2118956479, 1612420674, 1668583807, 1800004220, 1447372094,
  523904750, 1435821048, 923108080, 216161028, 1504871315, 306401572,
  2018281851, 1820959944, 2136819798, 359743094, 1354150250, 1843084537,
  1306570817, 244413420, 934220434, 672987810, 1686379655, 1301613820,
  1601294739, 484902984, 139978006, 503211273, 294184214, 176384212,
  281341425, 228223074, 147857043, 1893762099, 1896806882, 1947861263,
  1193650546, 273227984, 1236198663, 2116758626, 489389012, 593586330,
  275676551, 360187215, 267062626, 265012701, 719930310, 1621212876,
  2108097238, 2026501127, 1865626297, 894834024, 552005290, 1404522304,
  48964196, 5816381, 1889425288, 188942202, 509027654, 36125855,
  365326415, 790369079, 264348929, 513183458, 536647531, 13672163,
  313561074, 1730298077, 286900147, 1549759737, 1699573055, 776289160,
  2143346068, 1975249606, 1136476375, 262925046, 92778659, 1856406685,
  1884137923, 53392249, 1735424165, 1602280572
};

uint32_t dh_hash(const char *str)
{
  uint8_t c;
  uint32_t hash = 0;


  while( (c = *str++) )
    hash = (hash << 1) + dh_table[c];

  return hash;
}

/*
  The follwing hash algorithms appear to all come from
  Georgi Marinov
  Copyleft?
  http://www.sanmayce.com/Fastest_Hash/
 */

uint32_t alfalfa_hash(const char *str)
{
  size_t wrdlen = strlen(str);
  uint32_t hash = 7;


  for (uint32_t i = 0; i < (wrdlen & -2); i += 2) {
    hash = (17 + 9) * ((17 + 9) * hash + (str[i])) + (str[i+1]);
  }

  if (wrdlen & 1)
    hash = (17 + 9) * hash + (str[wrdlen-1]);

  return hash ^ (hash >> 16);
}

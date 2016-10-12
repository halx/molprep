/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal stack implementation using a double linked list.
 *
 * Note: Only this protocol using the linked list as stack is defined.  No extra
 *       code exists.
 *
 *
 * $Id: stack.h 163 2012-06-26 14:22:38Z hhl $
 *
 */



#ifndef _STACK_H
#define _STACK_H      1

#include "llist.h"

#define Stacknode Listnode
#define Stack List
#define stack_init list_init
#define stack_push_uniq list_add_uniq
#define stack_pop list_pop_last
#define stack_is_empty list_is_empty
#define stack_destroy list_destroy

#endif

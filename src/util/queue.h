/*
 * Copyright (C) 2012 Hannes Loeffler, STFC Daresbury, UK
 *
 *
 * A minimal queue/deque implementation using a double linked list.
 *
 * Note: Only this protocol using the linked list as queue is defined.  No extra
 *       code exists.
 *
 *
 * $Id: queue.h 163 2012-06-26 14:22:38Z hhl $
 *
 */



#ifndef _QUEUE_H
#define _QUEUE_H      1

#include "llist.h"

#define Queuenode Listnode
#define Queue List
#define queue_init list_init
#define queue_push_uniq list_add_uniq
#define queue_pop_front list_pop_first
#define queue_is_empty list_is_empty
#define queue_destroy list_destroy

#endif

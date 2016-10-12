#include <stdio.h>

#include "../util/llist.h"

int main(void)
{
  int nums[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  char *str[] = {"hello", "world", "bye, bye", "to", "you", "bla", "foo",
                 "bar"};

  char *data;

  List *list1 = NULL, *list2 = NULL;
  Listnode *node;



  list1 = list_init(list1);

  for (int i = 0; i < 8; i++) {
    list_add(list1, str[i]);
  }

  printf("nel = %i\n", list_num_el(list1) );

  FOREACH_LIST (node, list1) {
    printf("%s\n", (char *)list_get(node) );
  }

  printf("\n");

  list_destroy(list1);
  list1 = NULL;

  list2 = list_init(list2);

  for (int i = 0; i < 10; i++) {
    list_add(list2, &nums[i]);
  }

  printf("nel = %i\n", list_num_el(list2) );

  FOREACH_LIST (node, list2) {
    data = list_get(node);
    printf("%i\n", *data);
  }

  list_destroy(list2);
  list2 = NULL;
}

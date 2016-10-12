#include <stdio.h>

#include "../util/darray.h"

int main(void)
{
  int *data;
  int nums[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  char *str[] = {"hello", "world", "bye, bye", "to", "you", "bla", "foo"};

  DArray *array1 = NULL, *array2 = NULL;


  array1 = darray_init(array1, 3);

  for (int i = 0; i < 7; i++) {
    darray_add(array1, str[i]);
  }

  printf("nel = %i\n", darray_num_el(array1) );

  FOREACH_DARRAY (i, array1) {
    printf("%s\n", (char *)darray_get(array1, i) );
  }

  printf("\n");

  darray_destroy(array1);
  array1 = NULL;

  array2 = darray_init(array2, 3);

  for (int i = 0; i < 10; i++) {
    darray_add(array2, &nums[i]);
  }

  printf("nel = %i\n", darray_num_el(array2) );

  FOREACH_DARRAY (i, array2) {
    data = darray_get(array2, i);
    printf("%i\n", *data);
  }

  darray_destroy(array2);
  array2 = NULL;
}

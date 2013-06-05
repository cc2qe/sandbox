#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void plusEight(int *i) {

  *i = 8 + *i;

}


int main () {

  int *x = (int*)malloc(sizeof(int));
  *x = 20;

  plusEight(&x);
  printf("%d\n", *x);

  return 0;
}

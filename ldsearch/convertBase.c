#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int decToBase(int x,
	      int base)
{
  int i = 0;
  int r[64];
  char *buf = malloc(64);
  int x_b = 0;

  
  while(x != 0) {
    r[i] = x % base;
    x = x / base;
    ++i;
  }

  for (--i; i >= 0; --i) {
    sprintf(buf, "%s%d", buf, r[i]);
  }
  x_b = atoi(buf);

  return x_b;
}

int main()
{
  //int j;
  //for (j = 0; j <= 3; ++j) {
  //  printf("%d\n", j);
  //  printf("%d\n", decToBase(20,3));
  //}

  printf("%d\n", decToBase(5,3));
  printf("%d\n", decToBase(6,3));
  printf("%d\n", decToBase(7,3));
  printf("%d\n", decToBase(8,3));

  return 0;
}

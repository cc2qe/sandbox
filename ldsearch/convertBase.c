#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int decToBase(int x,
	      int base)
{
  int i = 0;
  int digit;
  int x_b = 0;
  
  while(x != 0) {
    digit = x % base;
    x = x / base;

    int tenPower = 1;
    int j;
    for (j = 0; j < i; ++j) {
      tenPower = tenPower * 10;
    }
    x_b += (digit * tenPower);

    ++i;
  }

  return x_b;
}

int main()
{
  int i;
  for (i = 0; i < 27; ++i) {
    printf("%03d\n", decToBase(i,3));
  }

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main() {
  char *str1 = "hi";
  char str2 = "a";

  char str3[100];

  int myLen = strlen(str3);

  printf("len is: %d\n", myLen);
  strcpy(str3, "TEXT");
  strcat(str3, str2);

  printf("%s\n", str3);

}

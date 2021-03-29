#include <stdio.h>
#include <string.h>

int main ()
{
  char str[80];
  strcpy (str,"strings  ");
  strcat (str,"have been ");
  strcat (str,"concatenated.");
  puts (str);
  return 0;
}

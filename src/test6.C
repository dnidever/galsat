#include <stdio.h>
#include <string.h>

int main ()
{
  char *stem = "galsat";
  char final[50];
  char str[7];
  int n=15;
  sprintf (str, "%d", n);
  //strcpy (stem,"galsat");
  strcat (final,stem);
  strcat (final,str);
  strcat (final,".out");
  puts (stem);
  puts (str);
  puts (final);

  //char buffer [50];
  //int n, a=5, b=3;
  //n=sprintf (buffer, "%d plus %d is %d", a, b, a+b);
  //printf ("[%s] is a %d chars string\n",buffer,n);

  return 0;
}

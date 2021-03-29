#include  <iostream>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
//#include  <unistd.h>                       // for getopt()
//#include  <fstream>                        // for file output
//#include  <stdio.h>
//#include  <stdlib.h>
//#include  <string.h>

int main()
{
  return 0;
}

//char * itoa(int n, char *buff, int radix)
// convert a positive integer n to char *buff
// for instant, this function work with radix <= 10;
// a little change to run with radix > 10
//{
char * buff;
int radix;
     int n = 5;
     int q;
     int r;
     int i = 0;
     char tmp[33];  // for radix = 2 and 32 bits computer

     q = int(n / radix);
//  while(q > 0) {
          q = int(n / radix);
          r = n % radix;
          n = q;
          tmp[i++] = 48 + r;
//   }
//while(q > 0);
     //int j;
     //for(j = 0; j < i; j++){
     //     buff[j] = tmp[i - j - 1];
     // }
     //buff[j] = NULL;
     //return buff;
//     return "A";
//}

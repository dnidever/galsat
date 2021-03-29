#include <iostream>
#include <string>
#include <cstdlib>
//using namespace std;

void main()
{
int i=12; // number to be converted
string s; // string to recive the number
char b[50]; // buffer of chars
int radix=10; // 2:bin, 8:octal, 10:dec, 16:hex
s=itoa(i,b,radix);
cout << "int i=" << s << endl;
cout << "buffer=" << b << endl;
 return 0;
}

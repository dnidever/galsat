#include  <iostream>
#include  <cstdlib>                        // for atoi() and atof()
#include  <fstream>
#include  <string>
using namespace std;

int main ()
{

  int nput = 1;
  int n = 10;
  char *filename = "galsat";
  //sprintf( filename, "%n", nput);
  const char *temp = ".out";
  strcat( filename, temp);

  //cout << "Hello World!" << endl;

  ofstream file ( "dum.txt", ios::out );
    
      if (file.is_open())
      {
	//file << "This is a line.\n";
	//file << "This is another line.\n";
        file << filename << endl;
        file << temp << endl;
      }
      file.close();

  return 1;
}


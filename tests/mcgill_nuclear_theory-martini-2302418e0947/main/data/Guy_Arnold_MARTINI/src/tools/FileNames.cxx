#include<string>
#include <sstream>
#include <iostream>
#include "FileNames.h"
using namespace std;
    
FileNames::FileNames() 
{
}

FileNames::~FileNames() 
{
}

string FileNames::GetFN(int i) const 
/*returns file name number i if specified in list.
  if not it generates a standard file name 
  "outputi.dat", where i is the number requested.*/
{
  if(i<=FileName.size()) 
    {
      return FileName.at(i-1);
    }
  else 
    {
      string istring;
      stringstream out;
      out << "output" << i << ".dat";
      istring = out.str();
      cout << "[FileNames speaking]: WARNING - file name not in list. Using standard file name: " << istring << endl;
      return istring;
    }
}

#ifndef FileNames_h //avoids multiple inclusions of the header file
#define FileNames_h
#include<string>
#include<vector> 

using namespace std;

class FileNames
{
 public:

  FileNames();
  ~FileNames();
  void SetFN(string x) {FileName.push_back(x);}
  string GetFN(int i) const;

 private:

  vector<string> FileName; 
};
#endif

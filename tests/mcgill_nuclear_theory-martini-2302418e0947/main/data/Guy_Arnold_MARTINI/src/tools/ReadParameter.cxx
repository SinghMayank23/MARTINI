#include "tools/ReadParameter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

ReadParameter::ReadParameter(string fn="param.inp")
{
    fname = fn;
}

ReadParameter::~ReadParameter()
{
    ifused.clear();
    linepos.clear();
    nameline.clear();
    varline.clear();
}

void ReadParameter::read()
{
  ifstream in;
  in.open(fname.c_str(),ios::in);
  if(!in) 
    {
      cerr << "Error. Unable to open file " << fname << endl;
      exit(1);
    }
  string templine;
  lineposition=0;
  while(getline(in,templine)) 
    {
      lineposition++;
      if(getParam(templine)) break;
    }
  in.close();
}

int ReadParameter::getParam(string templine)
{
  int icomment1=templine.find('#');
  int icomment2=templine.find('!');
  int icomment;

  if(icomment1*icomment2 > 0)
    icomment=min(icomment1,icomment2);
  else
    icomment=max(icomment1,icomment2);
  
  int ii=templine.find('=');
  
  if((icomment < 0 || icomment > ii) && ii >= 2) 
    {
      linepos.push_back(lineposition);
      
      nameline.push_back(templine.substr(0,ii));
      removeBlank(nameline.back());
      
      varline.push_back(templine.substr(ii+1));
      removeBlank(varline.back());
      
    } 
  else 
    {
      removeBlank(templine);
      if(templine == "END" || templine == "end" || templine == "End") 
	{
	  return 1;
	}
    }
  return 0;
}

void ReadParameter::removeBlank(string& line)
{
  const string delims(" \t,;");
  
  string::size_type begIdx, endIdx;
  begIdx = line.find_first_not_of(delims);
  line.erase(0,begIdx);
  
  // search end of the actural word
  endIdx = line.find_first_of(delims,begIdx);
  
  // end of word is end of line.
  if(endIdx == string::npos) endIdx = line.length();
  line.erase(endIdx);
}


//...Input integer variable from the config.
void ReadParameter::inputParam(string name,int& var)
{
  int nfound=0;
  for(int i=0; i< (int)nameline.size();i++) 
    {
      if(nameline[i] == name) 
	{
	  istringstream is(varline[i]);
	  if(is >> var) 
	    {
	      nfound++;
	      ifused.push_back(nfound);
	    }
	  else
	    {
	      cerr << "(::inputParam int) bad format " << nameline[i]
		   << " var = " << varline[i] << endl;
	    }
	}
    } 
}

//...Input floating variable from the config.
void ReadParameter::inputParam(string name,double& var)
{
  int nfound=0;
  for(int i=0; i<(int)nameline.size();i++) 
    {
      if(nameline[i] == name)
	{
	  istringstream is(varline[i]);
	  if(is >> var) 
	    {
	      nfound++;
	      ifused.push_back(nfound);
	    } 
	  else 
	    {
	      cerr << "(::inputParam double) bad format " << nameline[i]
		   << " var = " << varline[i] << endl;
	    }
	}
    }
}

//...Input string variable from the config.
void ReadParameter::inputParam(string name,string& var)
{
  int nfound=0;
  for(int i=0; i<(int)nameline.size();i++) 
    {
      if(nameline[i] == name) 
	{
	  var = varline[i];
	  nfound++;
	  ifused.push_back(nfound);
	}
    }
}

//...Check errors in the config file.
void ReadParameter::checkError()
{
  int nerr=0;
  for(int i=0;i<(int)ifused.size();i++)
    {
      if(ifused[i] != 1) 
	{
	  nerr++;
          if(ifused[i] == 0) 
	    {
	      cerr << "unknown " << nameline[i] << " at line "
		   << linepos[i] << endl;
	    } 
	  else 
	    {
	      cerr << "doubly defined [" << nameline[i]
		   << "] at line " << linepos[i] << endl;
	    }
	}
    }
  
  cout << "# " <<  ifused.size() << " parameter(s) set. " << nerr
       << " error(s) in input file" << endl;
  
  if(nerr > 0 ) exit(1);
}

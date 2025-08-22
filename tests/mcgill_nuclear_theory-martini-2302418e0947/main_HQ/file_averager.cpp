//I am finally making a file averaging program that accepts arguments and determines the range of 
//the files. -CFY 6/30/2011
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

//When calling this routine, it MUST have the arguments like so:
// ./file_averager name_tag_ min_index max_index n_columns
int main(int argc, char *argv[]){
  stringstream nameTagReader;
  nameTagReader << argv[1];
  string nameTag = nameTagReader.str();

  stringstream minReader;
  minReader << argv[2];
  string minString = minReader.str();
  int MIN = atoi(minString.c_str());

  stringstream maxReader;
  maxReader << argv[3];
  string maxString = maxReader.str();
  int MAX = atoi(maxString.c_str());

  stringstream nColumnsReader;
  nColumnsReader << argv[4];
  string nColumnsString = nColumnsReader.str();
  int nColumns = atoi(nColumnsString.c_str());

  int NFILES = MAX-MIN+1;
  cout << "Number of files to be read = " << NFILES << endl;
  double *x;
  double *y[nColumns];
  int NPOINTS = 0;
  bool findNPOINTS = true;

  for(int ifile = MIN; ifile <= MAX; ifile++){
    cout << ifile << endl;
    stringstream fileNumberReader;
    fileNumberReader << ifile;
    string fileNumber = fileNumberReader.str();
    string fileName = nameTag+fileNumber+".dat";

    ifstream file;
    file.open(fileName.c_str());
    if(!file){
      cout << fileName << " not found!" << endl;
      NFILES--;
    }
    else{
      if(findNPOINTS){
	ifstream justcounting(fileName.c_str());
	double ddummy;
	while(!justcounting.eof()){
	  for(int ic=0; ic<=nColumns; ic++){
	    justcounting >> ddummy;
	  }
	  if(!justcounting.eof()) NPOINTS++;
	}
	justcounting.close();

	x = new double[NPOINTS];
	for(int ic=0; ic<nColumns; ic++){
	  y[ic] = new double[NPOINTS];
	}

	findNPOINTS = false;
	cout << "NPOINTS = " << NPOINTS << endl;

	for(int ip=0; ip<NPOINTS; ip++){
	  x[ip] = 0.;
	  for(int ic=0; ic<nColumns; ic++){
	    y[ic][ip] = 0.;
	  }
	}
      }

      for(int ip=0; ip<NPOINTS; ip++){
	file >> x[ip];
	for(int ic=0; ic<nColumns; ic++){
	  double yValue;
	  file >> yValue;
	  y[ic][ip] += yValue;
	}
      }
    }
    file.close();
  }

  cout << "NFILES = " << NFILES << endl;

  for(int ip=0; ip<NPOINTS; ip++){
    for(int ic=0; ic<nColumns; ic++){
      y[ic][ip] /= (double)NFILES;
    }
  }

  string averagedFileName = nameTag+"averaged.dat";
  ofstream averagedFile;
  averagedFile.open(averagedFileName.c_str());

  for(int ip=0; ip<NPOINTS; ip++){
    averagedFile << x[ip] ;
    for(int ic=0; ic<nColumns; ic++){
      averagedFile << " " << y[ic][ip] ;
    }
    averagedFile << endl ;
  }

  averagedFile.close();

  return 0;
}

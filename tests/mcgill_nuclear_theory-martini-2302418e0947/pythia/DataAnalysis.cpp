#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#define nbins 100
#define ptmax 10.
#define size 10^5

using namespace std;

char *char_malloc(int n1);

int main() {

  string file_input = "./EventOutput/ALICE_pp.dat";
  fstream fin( file_input.c_str() ,ios::in); 

  int id[size];
  double pT[size];
  double y[size];
  double eta[size];

  int ival;
  double dval;

  int ik = 0;
  char *s;

  s = char_malloc(120);
  while ( strcmp(s, "EndOfEvent") != 0 ) {

    fin >> ival;
    fin >> dval;
    fin >> dval;
    fin >> dval;

    id[ik] = ival;
    pT[ik] = dval;
    y[ik] = dval;
    eta[ik] = dval;

    ik++;
  }

  cout << "ik = " << ik << endl;
}


char *char_malloc(int n1)
{
    char *char_ptr;

    /* pointer to the n1 array */
    char_ptr = (char *) malloc (sizeof(char)*n1);
}


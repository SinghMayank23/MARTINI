// Information.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// It has been modified from the PYTHIA event generator, which is
// Copyright (C) 2008 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains a class that keep track of generic event info.
// Information: contains information on the generation process and errors.

// Function definitions (not found in the header) for the Information class.

#include "Information.h"


//**************************************************************************

// Information class.
// This class contains a mixed bag of information on the event generation 
// activity, especially on the current subprocess properties.

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times the same error message will be repeated at most.
const int Information::TIMESTOPRINT = 1; 

//*********

// List (almost) all information currently set.

void Information::list(ostream& os) 
{


}

//*********
  
// Print a message the first few times. Insert in database.
 
void Information::errorMsg(string messageIn, string extraIn, ostream& os) {
   
  // Recover number of times message occured. Also inserts new string.
  int times = messages[messageIn];
  ++messages[messageIn];

  // Print message the first few times.
  if (times < TIMESTOPRINT) os << " MARTINI " << messageIn << " " 
    << extraIn << endl;

}

//*********

// Provide total number of errors/aborts/warnings experienced to date.

int Information::errorTotalNumber() {

  int nTot = 0;
  for ( map<string, int>::iterator messageEntry = messages.begin();
    messageEntry != messages.end(); ++messageEntry)
    nTot += messageEntry->second;
  return nTot;

}

//*********

// Print statistics on errors/aborts/warnings.

void Information::errorStatistics(ostream& os) {

  // Header.
  os << "\n *-------  MARTINI Error and Warning Messages Statistics  "
     << "-----------------------------------------------------------* \n"
     << " |                                                       "
     << "                                                          | \n"
     << " |  times   message                                      "
     << "                                                          | \n" 
     << " |                                                       "
     << "                                                          | \n";

  // Loop over all messages
  map<string, int>::iterator messageEntry = messages.begin();
  if (messageEntry == messages.end()) 
    os << " |      0   no errors or warnings to report              "
       << "                                                          | \n";
  while (messageEntry != messages.end()) {
    // Debug printout.
    string temp = messageEntry->first;
    int len = temp.length();
    temp.insert( len, max(0, 102 - len), ' ');
    os << " | " << setw(6) << messageEntry->second << "   " 
       << temp << " | \n";
    ++messageEntry;
  } 

  // Done. 
  os << " |                                                       "
     << "                                                          | \n"
     << " *-------  End MARTINI Error and Warning Messages Statistics"
     << "  -------------------------------------------------------* " 
     << endl;

}

//**************************************************************************

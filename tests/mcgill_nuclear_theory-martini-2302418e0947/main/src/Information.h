// Information.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// It has been modified from the PYTHIA event generator, which is
// Copyright (C) 2008 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains a class that keep track of generic event info.
// Information: contains information on the generation process and errors.

#ifndef Information_H
#define Information_H

#include "Pythia8/PythiaStdlib.h"
 
//**************************************************************************

// The Information class contains a mixed bag of information on the event
// generation activity, especially on the current subprocess properties,
// and on the number of errors encountered. This is used by the 
// generation machinery, but can also be read by the user.

using namespace Pythia8;

class Information {

public:

  // Constructor. 
  Information() {} 

  // Listing of most available information on current event.
  void   list(ostream& os = cout);
  
  // Reset to empty map of error messages.
  void errorReset() {messages.clear();}
  
  // Print a message the first few times. Insert in database.
  void errorMsg(string messageIn, string extraIn = " ",
    ostream& os = cout);

  // Provide total number of errors/aborts/warnings experienced to date.
  int  errorTotalNumber();

  // Print statistics on errors/aborts/warnings.
  void errorStatistics(ostream& os = cout);

private:

  // Friend classes allowed to set info.
  friend class Pythia;
  friend class ProcessLevel;
  friend class ProcessContainer;
  friend class PartonLevel;
  friend class MultipleInteractions;

  // Number of times the same error message is repeated.
  static const int TIMESTOPRINT;

  // Map for all error messages.
  map<string, int> messages;

};
 
//**************************************************************************

#endif // Information_H

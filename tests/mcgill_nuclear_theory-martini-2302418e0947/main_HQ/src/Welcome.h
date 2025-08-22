// Welcome.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the MARTINI welcome banner 

#ifndef Welcome_h //avoids multiple inclusions of the header file
#define Welcome_h

using namespace std;

class Welcome
{
 public:
  //constructor - display the welcome screen
  Welcome()
    {
      cout << "   ____________________________________________________ _____________________________ " << endl
	   << "  |                                                    |                             |" << endl
	   << "  |       __     __)                        \\__\\__/    | Main author:                |"  << endl 
	   << "  |      (, /|  /|              ,      ,     \\  O/     | Bjoern Schenke              |"  << endl 
	   << "  |        / | / |  _   __ _/_   __           \\ /      | McGill University, Montreal |"  << endl
	   << "  |     ) /  |/  |_(_(_/ (_(___(_/ (__(_       |       | Canada                      | "   << endl
	   << "  |    (_/   '                                _|_      |                             | " << endl
	   << "  |____________________________________________________|_____________________________|" << endl
	   << "  |                                                                                  |" << endl
	   << "  |     Main reference: B.Schenke, S. Jeon, C. Gale, Phys.Rev.C80, 054913 (2009)     |" << endl
	   << "  |__________________________________________________________________________________|" << endl
	   << " " << endl;       
    };
};  
#endif

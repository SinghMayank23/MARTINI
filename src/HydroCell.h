// HydroCell.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the HydroCell class which contains information on the soft background in on cell

#ifndef HydroCell_h
#define HydroCell_h

struct HydroCell
{
    double T;
    double vx;
    double vy;
    double vz;
    double QGPfrac;
};

#endif

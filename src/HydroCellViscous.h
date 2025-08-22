// HydroCell.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the HydroCell class which contains information on the soft background in on cell

#ifndef HydroCellViscous_h
#define HydroCellViscous_h

struct HydroCellViscous
{
    double Wtautau;
    double Wtaux;
    double Wtauy;
    double Wtaueta;
    double Wxx;
    double Wxy;
    double Wxeta;
    double Wyy;
    double Wyeta;
    double Wetaeta;
    double BulkPi;
    double Entropy;
    double cs2;
};

#endif

/**
 **        FILE: angles.h
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **              This header file contains the declarations for
 **              swiss.c and angles.c
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 **
 **  DISCLAIMER:   This program is free software; you can redistribute
 **              it and/or modify it under the terms of the GNU General
 **              Public License as published by the Free Software
 **              Foundation; either version 2 of the License, or
 **              (at your option) any later version.
 **                This program is distributed in the hope that it will
 **              be useful, but WITHOUT ANY WARRANTY; without even the
 **              implied warranty of MERCHANTABILITY or FITNESS FOR A
 **              PARTICULAR PURPOSE. See the GNU General Public License
 **              for more details.
 **                You should have received a copy of the GNU General
 **              Public License along with this program; if not, write
 **              to the Free Software Foundation, Inc., 675 Mass Ave,
 **              Cambridge, MA 02139, USA.
 **
 **      AUTHOR: Etienne Favey <eti@favey.lugs.ch>
 **              Zurich, Switzerland
 **
 **     VERSION: 4
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              03-12-1996 -- first release (Version 1.05)
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 **              07-05-1997 -- added admsstr(), ellsupport, det_ell()
 **              14-01-1998 -- added struct frameconv, select_ellipsoid(),
 **                            select_frameconv()
 **              20-01-1998 -- split swiss.c into swiss.c and trafo.c
 **              13-05-1998 -- parameter for ch_magnetic_deviation changed
 **
 **/
#ifndef ANGLES
#define ANGLES


/* ----------------------------------------------------------------
 * for angles.c
 */

struct hms {
    int deg;
    int min;
    double sec;
  } ;

struct Angle {
  double rad;
  double grad;
  double deg;
  struct hms dms;
} ;

struct Angle afromrad (double angrad);
struct Angle afromdeg (double angdeg);
struct Angle afromgrad(double anggrad);

struct Angle afromdmsverb(int degrees, int minutes, double seconds);
struct Angle afromdms(struct hms mydms);
struct Angle afromstr(char *s);

char * admsstr (const struct Angle a, const int mode);
/*
 * mode: 0 -> 46 30 25.643346789
 *       1 -> 46d30m25.643346789
 *       2 -> 46\t30\t25.643346789
 */
#endif


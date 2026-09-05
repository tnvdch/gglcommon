/**
 **        FILE: swiss.h
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
#ifndef SWISS
#define SWISS

/* ----------------------------------------------------------------
 * for swiss.c
 */

int ch_elltosph(double elllat, double elllon, double *sphlat, double *sphlon);
/*
 *  conform projection from ellipsoid to sphere
 *  parameters have to be passed and are returned in rad
 */

int ch_sphtoell(double sphlat,double sphlon, double *elllat, double *elllon);
/*
 *  conform reverse projection from sphere to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

int ch_scale_etom(double elllat, double sphlat, double pselat,double *scale);
/*
 *  scale factor of projection from ellipsiod to plane
 *  is only dependent on latitude
 */

int ch_sphtopse(double sphlat, double sphlon, double *pselat, double *pselon);
/* 
 *  projection from sphere to pseudo equator system
 *  parameters have to be passed and are returned in rad
 */

int ch_psetosph(double pselat, double pselon, double *sphlat, double *sphlon);
/* 
 *  reverse projection from pseudo to sphere
 *  parameters have to be passed and are returned in rad
 */

int ch_merconv(double sphlat, double sphlon, double *conv);
/*
 *  meridian convergence for swiss projection
 */

int ch_psetomil(double pselat, double pselon, double *milx, double *mily);
/* 
 *  mercator projection pseudo to plane
 *  parameters have to be passed and are returned in rad
 */

int ch_miltopse(double milx, double mily, double *pselat, double *pselon);
/*
 *  reverse mercator projection plane to pseudo
 *  parameters have to be passed and are returned in rad
 */

int ch_elltomil(double elllat, double elllon, double *milx, double *mily);
/*
 *  conform projection from ellipsoid to plane (mil)
 *  parameters have to be passed and are returned in rad
 */


int ch_miltoell(double milx, double mily, double *elllat, double *elllon);
/*
 *  reverse mercator projection plane to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

int ch_e2m(double elllat, double elllon, 
	   double *milx, double *mily, 
	   double *scale, double *merconv);
/*
 *  conform projection from ellipsoid to plane (mil)
 *  parameters have to be passed and are returned in rad
 */

int ch_m2e(double milx, double mily,
	   double *elllat, double *elllon,
	   double *scale, double *merconv);
/* 
 *  reverse mercator projection plane to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

double ch_magnetic_deviation (
  const double milx,
  const double mily,
  const double year);

#endif


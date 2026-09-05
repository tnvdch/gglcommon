/**
 **        FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/swiss.c,v $
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 ** $Id: swiss.c,v 1.5 1999/01/26 17:19:19 favey Exp $
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
 **     CHANGES: dd-mm-yyyy -- description
 **              27-08-1996 -- first release (Version 1.01)
 **              03-12-1996 -- included xyz_to_blh, blh_to_xyz
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 ** $Log: swiss.c,v $
 ** Revision 1.5  1999/01/26 17:19:19  favey
 ** moved out of swiss.h
 **
 ** Revision 1.4  1998/10/30 14:22:07  favey
 ** changed config.h to global.h
 **
 ** Revision 1.3  1998/05/13 13:54:52  cocard
 ** *** empty log message ***
 **
 ** Revision 1.2  1998/05/06 09:42:35  cocard
 ** old fashioned function headers replaced (MC)
 **
 ** Revision 1.1  1998/04/22 12:35:09  favey
 ** Initial revision
 **
 ** Revision 0.6  1998/03/03 10:37:45  favey
 ** added ch_magnetic_deviation
 **
 ** Revision 0.5  1998/01/20 09:59:00  favey
 ** split swiss.c into swiss.c and trafo.c
 **
 ** Revision 0.4  1998/01/15 17:22:24  favey
 ** new functions select_[<ellipsoid><frameconv><formtype>]
 **
 ** Revision 0.3  1997/12/24 12:12:27  favey
 ** followup of V. 0.29
 **
 **
 **/
/*
 *  the standard C headers
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

/*
 *  the GGL Common Library headers
 */
#include "global.h"
#include "errproc.h"
#include "angles.h"
#include "swiss.h"

/*
 * Swiss Datum only! used by swiss.c
 *
 * Bessel Ellisoid as used for the swiss projection
 */
#define BESSEL_A 6377397.15500
#define BESSEL_B 6356078.96282
#define BESSEL_SQ_E1 0.006674372230614
#define BESSEL_SQ_E2 0.006719218797971

/*
 * Reference Point (old observatory in Bern) on ellipsoid
 */
#define BERNOBS_DEG_LAT 46.9524055556
#define BERNOBS_DEG_LON  7.4395833333
#define BERNOBS_RAD_LAT  0.819474068677
#define BERNOBS_RAD_LON  0.129845224143

/*
 * Curvature Radii in Reference Point
 */
#define BERNOBS_W 0.998216186318
#define BERNOBS_RM 6368853.81749
#define BERNOBS_RN 6388793.57238
#define BERNOBS_R0 6378815.90365

/*
 * Constants for conform projection Ellipsoid-Sphere after Gauss
 */
#define ELLSPH_DEG_B0 46.907731457905
#define ELLSPH_B0 0.81869435858
#define ELLSPH_AL 1.0007291384304
#define ELLSPH_K  1.0030714396117

/*
 * Macros for Ellipsoid calculations
 */
#define ELLCALC_W(sqre,lat) (sqrt(1 - sqre * sin(lat) * sin(lat)))



/**
 ** Offset for military coordinates
 **/
#define OFF_MIL_X 200000
#define OFF_MIL_Y 600000

/**
 ** Criterion to abort iteration
 **/
#define TOLERANCE 10e-10

/**
 ** conform projection from ellipsoid to sphere
 ** parameters have to be passed and are returned in rad
 **/
int ch_elltosph(double elllat, double elllon, double *sphlat, double *sphlon)
{
  double e1, pex;

  e1 = sqrt(BESSEL_SQ_E1);

  pex = ELLSPH_AL * log(tan(M_PI/4.0 + elllat/2.0))
    + ELLSPH_AL * e1 / 2 * log((1-e1*sin(elllat))/(1+e1*sin(elllat)))
    + log(ELLSPH_K);

  *sphlat = 2*(atan(exp(pex)) - M_PI/4);
  *sphlon = ELLSPH_AL * (elllon - BERNOBS_RAD_LON);

  return EXIT_SUCCESS;
}


/**
 ** conform reverse projection from sphere to ellipsoid
 ** parameters have to be passed and are returned in rad
 **/
int ch_sphtoell(double sphlat,double sphlon, double *elllat, double *elllon)
{
  double e1, brackets;
  double latiter, olditer;

  e1 = sqrt(BESSEL_SQ_E1);
  latiter = sphlat;

  do
    {
      brackets = pow((tan(M_PI/4.0 + sphlat/2.0) / ELLSPH_K),(1.0/ELLSPH_AL))
	* pow(((1.0+e1*sin(olditer=latiter))/(1.0-e1*sin(latiter))),(e1/2.0));

      latiter = 2.0*atan(brackets) - M_PI/2.0;
    } while (fabs(latiter-olditer) > TOLERANCE);
  *elllat = latiter;
  *elllon = BERNOBS_RAD_LON + sphlon/ELLSPH_AL;

  return EXIT_SUCCESS;
}

/**
 ** scale factor of projection from ellipsiod to plane
 ** is only dependent on latitude
 **/
int ch_scale_etom(double elllat, double sphlat, double pselat,double *scale)
{
  double radius_N;

  radius_N = BESSEL_A / ELLCALC_W(BESSEL_SQ_E1,elllat);

  *scale = ELLSPH_AL * BERNOBS_R0 * cos(sphlat) /
    ( radius_N * cos(elllat) * cos(pselat) );

  return EXIT_SUCCESS;
}

/**
 ** projection from sphere to pseudo equator system
 ** parameters have to be passed and are returned in rad
 **/
int ch_sphtopse(double sphlat, double sphlon, double *pselat, double *pselon)
{
  *pselon = atan2(sin(sphlon), 
		  (sin(ELLSPH_B0)*tan(sphlat) + cos(ELLSPH_B0)*cos(sphlon)) );

  *pselat = asin(cos(ELLSPH_B0)*sin(sphlat) 
		 - sin(ELLSPH_B0)*cos(sphlat)*cos(sphlon) );

  return EXIT_SUCCESS;
}

/**
 ** reverse projection from pseudo to sphere
 ** parameters have to be passed and are returned in rad
 **/
int ch_psetosph(double pselat, double pselon, double *sphlat, double *sphlon)
{
  *sphlat = asin(cos(ELLSPH_B0)*sin(pselat) 
		 + sin(ELLSPH_B0)*cos(pselat)*cos(pselon) );

  *sphlon = atan2(sin(pselon),
		  (cos(ELLSPH_B0)*cos(pselon) - sin(ELLSPH_B0)*tan(pselat)) );

  return EXIT_SUCCESS;
}

/**
 ** meridian convergence for swiss projection
 **/
int ch_merconv(double sphlat, double sphlon, double *conv)
{
  *conv = atan2(sin(sphlon),
		(cos(sphlat)/tan(ELLSPH_B0) + sin(sphlat)*cos(sphlon)) );

  return EXIT_SUCCESS;
}

/**
 ** mercator projection pseudo to plane
 ** parameters have to be passed and are returned in rad
 **/
int ch_psetomil(double pselat, double pselon, double *milx, double *mily)
{
  *milx = BERNOBS_R0 * log(tan(M_PI/4 + pselat/2)) + OFF_MIL_X;

  *mily = BERNOBS_R0 * pselon + OFF_MIL_Y;

  return EXIT_SUCCESS;
}

/**
 ** reverse mercator projection plane to pseudo
 ** parameters have to be passed and are returned in rad
 **/
int ch_miltopse(double milx, double mily, double *pselat, double *pselon)
{
  *pselat = 2 * atan(exp((milx-OFF_MIL_X)/BERNOBS_R0)) - M_PI/2;

  *pselon = (mily-OFF_MIL_Y)/BERNOBS_R0;

  return EXIT_SUCCESS;
}

/**
 ** conform projection from ellipsoid to plane (mil)
 ** parameters have to be passed and are returned in rad
 **/
int ch_elltomil(double elllat, double elllon, double *milx, double *mily)
{
  double sphlat, sphlon, pselat, pselon;

  ch_elltosph(elllat, elllon, &sphlat, &sphlon);
  ch_sphtopse(sphlat, sphlon, &pselat, &pselon);
  ch_psetomil(pselat, pselon, &(*milx), &(*mily));

  return EXIT_SUCCESS;
}

/**
 ** reverse mercator projection plane to ellipsoid
 ** parameters have to be passed and are returned in rad
 **/
int ch_miltoell(double milx, double mily, double *elllat, double *elllon)
{
  double sphlat, sphlon, pselat, pselon;

  ch_miltopse(milx, mily, &pselat, &pselon);
  ch_psetosph(pselat, pselon, &sphlat, &sphlon);
  ch_sphtoell(sphlat, sphlon, &(*elllat), &(*elllon));

  return EXIT_SUCCESS;
}

/**
 ** conform projection from ellipsoid to plane (mil)
 ** parameters have to be passed and are returned in rad
 **/
int ch_e2m(double elllat, double elllon, 
	   double *milx,  double *mily, 
	   double *scale, double *merconv)
{
  double sphlat, sphlon, pselat, pselon;

  ch_elltosph(elllat, elllon, &sphlat, &sphlon);
  ch_sphtopse(sphlat, sphlon, &pselat, &pselon);
  ch_psetomil(pselat, pselon, &(*milx), &(*mily));
  ch_merconv(sphlat, sphlon, &(*merconv));
  ch_scale_etom(elllat, sphlat, pselat, &(*scale));

  return EXIT_SUCCESS;
}

/**
 ** reverse mercator projection plane to ellipsoid
 ** parameters have to be passed and are returned in rad
 **/
int ch_m2e(double milx, double mily,
	   double *elllat, double *elllon,
	   double *scale, double *merconv)
{
  double sphlat, sphlon, pselat, pselon;

  ch_miltopse(milx, mily, &pselat, &pselon);
  ch_psetosph(pselat, pselon, &sphlat, &sphlon);
  ch_sphtoell(sphlat, sphlon, &(*elllat), &(*elllon));
  ch_merconv(sphlat, sphlon, &(*merconv));
  ch_scale_etom(*elllat, sphlat, pselat, &(*scale));

  return EXIT_SUCCESS;
}


double ch_magnetic_deviation (
  const double milx,
  const double mily,
  const double year)
{
  double dx, dy, dt;
  double d, g, f;

  dx = milx/1000.0 - 200.0;
  dy = mily/1000.0 - 600.0;
  dt = year - 1978;

  d = -2.6124 - 9.289e-4*dx + 4.9061e-3*dy;
  g = 0.16527*(1 + 4.6682e-4*dx - 6.36624e-4*dy);
  f = (dt + 1.1662e-2*(dt*dt))/(1 + 5.1943e-2*dt);
  d = d + g*f;
  return (d * M_PI / 180);
}


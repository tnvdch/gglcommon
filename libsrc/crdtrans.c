/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/crdtrans.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations on 3x3 matrix and 3x1 vectors
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: crdtrans.c,v 1.12 2001/07/03 08:28:49 cocard Exp favey $
 *
 *  DISCLAIMER:   This program is free software; you can redistribute
 *              it and/or modify it under the terms of the GNU General
 *              Public License as published by the Free Software
 *              Foundation; either version 2 of the License, or
 *              (at your option) any later version.
 *                This program is distributed in the hope that it will
 *              be useful, but WITHOUT ANY WARRANTY; without even the
 *              implied warranty of MERCHANTABILITY or FITNESS FOR A
 *              PARTICULAR PURPOSE. See the GNU General Public License
 *              for more details.
 *                You should have received a copy of the GNU General
 *              Public License along with this program; if not, write
 *              to the Free Software Foundation, Inc., 675 Mass Ave,
 *              Cambridge, MA 02139, USA.
 *
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *            Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/crdtrans.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.12  2001/07/03 08:28:49  cocard
 * *** empty log message ***
 *
 * Revision 1.11  2001/06/08 12:32:23  cocard
 * *** empty log message ***
 *
 * Revision 1.10  1999/10/07 14:52:17  cocard
 * *** empty log message ***
 *
 * Revision 1.9  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.8  1999/02/10 18:03:59  favey
 * changes to make commonlib compile on Solaris
 *
 * Revision 1.7  1999/02/10 17:27:52  favey
 * changed PI to M_PI for portability
 *
 * Revision 1.6  1999/01/26 17:19:19  favey
 * moved out of swiss.h
 *
 * Revision 1.5  1998/08/27 18:58:21  cocard
 * bug fixed in xyz_to_azel
 *
 * Revision 1.4  1998/08/20 08:36:49  cocard
 * center_angle numerically stabilized
 *
 * Revision 1.3  1998/06/17 15:53:40  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1998/06/15 15:38:15  cocard
 * function center_angle added
 *
 * Revision 1.1  1998/04/22 12:35:56  favey
 * Initial revision
 *
 * Revision 0.3  1997/12/24 12:08:55  favey
 * followup of V. 0.29
 *
 *
 *
 */

/*
 *  standard C libraries
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
 *  GGL Common Libs
 */
#include "global.h"
#include "trafo.h"
#include "crdtrans.h"
#include "mat3.h"

double mod2pi (double x)
{
  while (x>=2*M_PI) x -= 2*M_PI;
  while (x<0) x+= 2*M_PI;
  return x;
}

double ellipse_circumference (const double a,
			      const double ex)
{
  double x, dx, fac;
  int i = 0;
  fac = x = 1.0;
  do
    {
      i += 2;
      fac *= ex*(i-1)/(double) i;
      dx = fac*fac/(i-1);
      x -= dx;
    } while (dx>=1.0e-17);
  return 2*M_PI*a*x;
}

double center_angle( const double lat1,
		     const double lon1,
		     const double lat2,
		     const double lon2)
{
  double x;
  x = sin(lat1)*sin(lat2);
  x += cos(lat1)*cos(lat2)*cos(lon2-lon1);
  if (x> 1) return 0;
  if (x<-1) return M_PI;
  return acos(x);
}

double ellipse_approxdist (const double lat1,
			   const double lon1,
			   const double lat2,
			   const double lon2,
			   const struct ellipsoid ell)
{
  double w2, rn2, rm2, dlat, dlon, b3, r;

  b3 = (lat1+lat2)/2.0;
  dlat = lat2 - lat1;
  dlon = lon2 - lon1;
  dlon = mod2pi (dlon);
  if (dlon>M_PI) dlon = 2*M_PI - dlon;
  w2 = sin (b3);
  w2 = w2*w2;
  w2 = 1 - ell.e2*w2;

  rn2 = ell.a*ell.a/w2;

  w2 = w2*w2*w2;
  rm2 = ell.a*(1-ell.e2);
  rm2 = rm2*rm2/w2;

  r = cos (b3);
  r = r*r*rn2*dlon*dlon;
  r += rm2*dlat*dlat;
  return sqrt(r);
}

int xyz_to_azel (double xyz1[3],
		 double xyz2[3],
		 struct ellipsoid ell,
		 double *dist,
		 double *elev,
		 double *azi)
{
  double dxyz[3];
  double b, l, h;

  xyz_to_llh (ell, xyz1, &b, &l, &h);
  v3_sub (xyz2, xyz1, dxyz);
  *dist = v3_length (dxyz);
  rotate (dxyz, 3, l, dxyz);
  rotate (dxyz, 2, M_PI/2-b, dxyz);
  *elev = M_PI/2 - acos(dxyz[2]/ *dist);
  *azi = atan2(dxyz[1], -dxyz[0]);
  if (*azi < 0) (*azi) += 2*PI; 
  return 0;
}

int top_to_xyz (double lat,
		double lon,
		double dtop[3],
		double dxyz[3])
{
  rotate (dtop, 1, lat-M_PI/2, dxyz);
  rotate (dxyz, 3,-lon-M_PI/2, dxyz);
  return 0;
}

int top_to_xyz2 (double xyz[3],
		 double dtop[3],
		 double dxyz[3],
		 struct ellipsoid ell)
{
  double b,l,h;

  xyz_to_llh (ell, xyz, &b, &l, &h);
  top_to_xyz (b, l, dtop, dxyz);
  return 0;
}

int azel_to_xyz (double xyz[3],
		 double azi,
		 double ele,
		 double dxyz[3],
		 struct ellipsoid ell)
{
  double dtop[3];
  
  dtop[0] = cos(ele) * sin(azi);
  dtop[1] = cos(ele) * cos(azi);
  dtop[2] = sin(ele);
  top_to_xyz2 (xyz, dtop, dxyz, ell);
  return 0;
}

int xyz_to_top (double lat,
		double lon,
		double dxyz[3],
		double dtop[3])
{
  rotate (dxyz, 3, lon+M_PI/2, dtop);
  rotate (dtop, 1,-lat+M_PI/2, dtop);
  return 0;
}

int xyz_to_top2 (double xyz[3],
		 double dxyz[3],
		 double dtop[3],
		 struct ellipsoid ell)
{
  double b,l,h;

  xyz_to_llh (ell, xyz, &b, &l, &h);
  return xyz_to_top (b, l, dxyz, dtop);
}

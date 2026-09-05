/*
 *      FILE: crdtrans
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations on 3x3 matrix and 3x1 vectors
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *
 *
 */

#ifndef GGL_CRDTRANS
#define GGL_CRDTRANS


double mod2pi (double x)
;

double center_angle( const double lat1,
		     const double lon1,
		     const double lat2,
		     const double lon2)
;

double ellipse_circumference (const double a,
			      const double ex)
;

double ellipse_approxdist (const double lat1,
			   const double lon1,
			   const double lat2,
			   const double lon2,
			   const struct ellipsoid ell)
;

int xyz_to_azel (double xyz1[3],
		 double xyz2[3],
		 struct ellipsoid ell,
		 double *dist,
		 double *elev,
		 double *azi)
;

int top_to_xyz (double lat,
		double lon,
		double dtop[3],
		double dxyz[3])
;

int top_to_xyz2 (double xyz[3],
		 double dtop[3],
		 double dxyz[3],
		 struct ellipsoid ell)
;

int azel_to_xyz (double xyz[3],
		 double azi,
		 double ele,
		 double dxyz[3],
		 struct ellipsoid ell)
;

int xyz_to_top (double lat,
		double lon,
		double dxyz[3],
		double dtop[3])
;

int xyz_to_top2 (double xyz[3],
		 double dxyz[3],
		 double dtop[3],
		 struct ellipsoid ell)
;

#endif

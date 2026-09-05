/*
 *      FILE: intersect
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for geometrical intersections
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
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            15-05-1998  --  created
 *
 */

#ifndef INTERSECT
#define INTERSECT

/*
 * calculates the coordinates of a point on a line
 */
void line_coord
(const double x[3],           /* input:  coordinates of point */
 const double e[3],           /* input:  unary direction vector */
 const double dist,           /* input:  distance to new point */
 double xnew[3])              /* output: coordinates of new point */
;

/*
 * calculates the intersection of a line and an ellipsoid 
 * defined by its three axis
 * return value = #count of intersections (0,1, or 2)
 */
int intersection_line_ell3
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double axis[3],        /* input:  all 3 major axis  */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line and a ellipoid of rotation
 * return value = #count of intersections (0,1, or 2)
 * derivated from intersection_line_ell3
*/
int intersection_line_ell
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double major_axis,     /* input:  major axis */ 
 const double ex2,            /* input:  square value of excentricity */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line and a sphere
 * return value = #count of intersections (0,1, or 2)
 * derivated from intersection_line_ell3
 */
int intersection_line_sphere
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double radius,         /* input:  radius of sphere */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line with a surface with
 * constant longitude (plane surface)
 * return value = #count of intersections (-1, 0 or 1)
 */
int intersection_line_longitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double longitude,      /* input:  longitude */
 double *dist)                /* output: distance to intersection */
;

/*
 * calculates the intersection of a line with a surface with
 * constant elliposidal latitude (conical surface).
 * return value = #count of intersections (-1, 0, 1, or 2)
 */
int intersection_line_latitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double latitude,       /* input:  latitude */
 const double major_axis,     /* input:  major axis */ 
 const double ex2,            /* input:  excentricity of ellipsoid squared */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * 2D routine :
 * This function calculates the intersection of two circles given
 * by the coordinates  of the center and the radius. The resulting
 * intersection points are returned in res1 and res2. The return
 * value of the function specifies the number of valid inter-
 * sections (-1, 0, 1, 2) where -1 means circles are identical /
 * 0 = no intersection / 1 = one tangential point / 2 = 2 inter-
 * section points.
 */
int intersection_circles (
  const double center1[2],  /* input:  coordinates of center of the first circle */
  const double radius1,     /* input:  radius of the first circle */
  const double center2[2],  /* input:  coordinates of center of the second circle */
  const double radius2,     /* input:  radius of the second circle */
  double res1[2],           /* output: coordinates of intersection 1 */
  double res2[2])           /* output: coordinates of intersection 2 */
;

int intersection_line_segment(
  double origin[2], 
  double azi, 
  double p1[2], 
  double p2[2], 
  double ipoint[2])
;

#endif

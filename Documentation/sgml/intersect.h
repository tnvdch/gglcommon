/*
 *      FILE: ggl_time
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for geometrical intersections
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Cocard Marc     (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            15-05-1998  --  created
 *
 */
#ifndef INTERSECT
#define INTERSECT
/*
  calculates the coordinates of a point on a line
 */
void line_coord
(const double x[3],           /* input:  coordinates of point */
 const double e[3],           /* input:  unary direction vector */
 const double dist,           /* input:  distance to new point */
 double xnew[3])              /* output: coordinates of new point */
;

/*
  calculates the intersection of a line and a sphere
  return value = #count of intersections (0,1, or 2)
*/
int intersection_line_sphere
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double radius,         /* input:  radius of sphere */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
  calculates the intersection of a line with a surface with
  constant longitude (plane surface)
  return value = #count of intersections (-1, 0 or 1)
 */
int intersection_line_longitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double longitude,      /* input:  longitude */
 double *dist)                /* output: distance to intersection */
;

/*
  calculates the intersection of a line with a surface with
  constant latitude (conical surface)
  return value = #count of intersections (-1, 0, 1, or 2)
 */
int intersection_line_latitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double latitude,       /* input:  latitude */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;


int arc_intersect (
  const double p1[3],
  const double p2[3],
  double na[2],
  double nb[2])
;




#endif

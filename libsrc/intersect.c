/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/intersect.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for geometrical intersections
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: intersect.c,v 1.13 2001/07/27 17:18:15 cocard Exp favey $
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
 *   AUTHORS: Cocard Marc   <cocard@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            02-05-1998  --  created
 *
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/intersect.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.13  2001/07/27 17:18:15  cocard
 * *** empty log message ***
 *
 * Revision 1.12  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.11  1998/09/03 11:39:09  cocard
 * *** empty log message ***
 *
 * Revision 1.10  1998/09/03 09:39:35  cocard
 * intersection_line_latitude (bugfix)
 *
 * Revision 1.9  1998/09/03 09:34:39  cocard
 * *** empty log message ***
 *
 * Revision 1.8  1998/05/29 09:05:18  cocard
 * cosmetic
 *
 * Revision 1.7  1998/05/27 09:28:15  favey
 * removed include of common.h, not needed
 *
 * Revision 1.6  1998/05/20 09:04:58  cocard
 * new function intersect_circles added (MC)
 *
 * Revision 1.5  1998/05/15 19:49:37  cocard
 * *** empty log message ***
 *
 * Revision 1.4  1998/05/15 19:47:44  cocard
 * first check in 15/05/98 (MC)
 *
 */

/*
 *  Standard Libraries
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "intersect.h"


/*
 * 3D: calculates the coordinates of a point on a line 
 */
void line_coord
(const double x[3],           /* input:  coordinates of point */
 const double e[3],           /* input:  unary direction vector */
 const double dist,           /* input:  distance to new point */
 double xnew[3])              /* output: coordinates of new point */
{
  int i;
  for (i=0;i<3;i++) xnew[i]=x[i]+dist*e[i];
}


/*
 * 3D: calculates the intersection of a line and an ellipsoid
 * defined by its three axis
 * return value = #count of intersections (0,1, or 2)
 */
int intersection_line_ell3
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double axis[3],        /* input:  all 3 major axis  */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
{
  int i;
  double a, b, c, det;
  double xr, er;
  a=0;  b=0; c=-1;
  for (i=0; i<3; i++)
    {
      xr=x[i]/axis[i];
      er=e[i]/axis[i];
      a += er*er;
      b += er*xr;
      c += xr*xr;
    }
  
  det = b*b-a*c;
  if (det<0)  return 0;           /* no intersection */ 
  det=sqrt(det);
  *dist1=(-b-det)/a;
  *dist2=(-b+det)/a;
  if (*dist1==*dist2) return 1;   /* same point */
  return 2;
}

/*
 * 3D: calculates the intersection of a line and a ellipoid of rotation 
 * return value = #count of intersections (0,1, or 2)
 */
int intersection_line_ell
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double major_axis,     /* input:  major axis */ 
 const double ex2,            /* input:  square value of excentricity */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
{
  double axis[3];
  axis[0]=axis[1]=major_axis;
  axis[2]=sqrt(1-ex2) * major_axis;
  return intersection_line_ell3(x, e, axis, dist1, dist2);
}

/*
 * 3D: calculates the intersection of a line and a sphere
 * return value = #count of intersections (0,1, or 2)
 */
int intersection_line_sphere
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double radius,         /* input:  radius of sphere */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
{
  double axis[3];
  axis[0]=axis[1]=axis[2]=radius;
  return intersection_line_ell3(x, e, axis, dist1, dist2);
}

/*
 * 3D: calculates the intersection of a line with a surface with 
 * constant longitude (plane surface) 
 * return value = #count of intersections (-1, 0 or 1)
 */
int intersection_line_longitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double longitude,      /* input:  longitude */
 double *dist)                /* output: distance to intersection */
{
  double co, si, denominator;

  co=cos(longitude);
  si=sin(longitude);
  if ((denominator=e[1]*co-e[0]*si) == 0)
    { 
      if (x[1]*co-x[0]*si == 0) return -1;  /* all points on line are ok */
      return 0;                             /* no intersection */
    }
  *dist=-(x[1]*co-x[0]*si)/denominator;
  return 1;
}

/*
 * 3D: calculates the intersection of a line with a surface with
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
{
  const double zero = 1.0e-12;

  double si, sin2, cos2;      
  double a,b,c;               
  double det;
  double eps, x2, dummy;
  int intersection;

  si=sin(latitude);    
  sin2 = si*si;  
  cos2 = 1 - sin2;

  eps = major_axis/sqrt(1-ex2*sin2);
  eps *= ex2*si;
  x2 = x[2] + eps;

  a = cos2*e[2]*e[2] - sin2*(e[0]*e[0]+e[1]*e[1]);
  b = cos2* x2 *e[2] - sin2*(x[0]*e[0]+x[1]*e[1]);
  c = cos2* x2 * x2  - sin2*(x[0]*x[0]+x[1]*x[1]);
  det=b*b - a*c;

  if (det<0) return 0;                  /* no intersections */

  if (fabs(a)<zero)                     /* degrades to a  linear equation */
    {
      if (fabs(b)<zero)
        {
          if (fabs(c)<zero) return -1;  /* all points on line are ok */
          return 0;                     /* else no intersection */
        }
      *dist1 = -c/2/b;                  /* case where b not zero */
      if ((x[2]+(*dist1)*e[2]) * si >= 0)  return 1;
      return 0;
    }

  det=sqrt(det);
  *dist1=(-b+det)/a;
  *dist2=(-b-det)/a;
 
  intersection=0;
  if ((x[2]+(*dist1)*e[2])*si >= 0) intersection++; else *dist1=*dist2;
  if ((x[2]+(*dist2)*e[2])*si >= 0) intersection++; else *dist2=*dist1;

  if (*dist2<*dist1)
  {
    dummy  = *dist2; 
    *dist2 = *dist1; 
    *dist1 = dummy;
  }
  return intersection;
}

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
{
  double dcenter[2], dx[2], e[2];
  double dist, dist2;
  double dummy;
  int i;

#ifdef WITH_ERRORS
    if (radius1==0 || radius2==0)
      {
	fprintf (stderr, "GGL Common Library --- intersect\n");
	fprintf (stderr, "intersection_circles: radius of circle = 0\n");
      }
#endif

  /* calculation of distance between centers */
  dist2=0;
  for (i=0;i<2;i++) 
    {
      dcenter[i]=center2[i]-center1[i];
      dist2 += dcenter[i]*dcenter[i];
    }
  dist = sqrt(dist2);
  if (dist==0)
    {
      if (radius1==radius2) return -1;  /* circles identical */
      return 0;                         /* no intersection */
    }
  /* cosine equation */
  dummy=radius1*radius1;
  dx[0] = (dist2 + dummy - radius2*radius2)/(2*dist);
  dx[1] = dummy - dx[0]*dx[0];

  if (dx[1]<0) return 0;                /* no intersection */
  dx[1]=sqrt(dx[1]);

  e[0]=dcenter[0]/dist;
  e[1]=dcenter[1]/dist;

  res1[0] = res2[0] = center1[0] + dx[0]*e[0];
  res1[1] = res2[1] = center1[1] + dx[0]*e[1];

  res1[0] -= dx[1]*e[1];  res1[1] += dx[1]*e[0]; /* first solution */
  res2[0] += dx[1]*e[1];  res2[1] -= dx[1]*e[0]; /* second solution */

  if (dx[1]==0) return 1;           /* 1 tangential point */
  return 2;                         /* 2 intersection points */
}

/* 
 * 2D: calculates the intersection of a line given by
 * its origin and the direction (azimuth)
 * with a segment defined by 2 points
 * returns 0 if intersecting -1 if no intersection 
 * if ipoint != NULL the intersection point is stored there
 * -1 if no intersection 
 */

int intersection_line_segment(double origin[2], double azi, double p1[2], double p2[2], double ipoint[2])
{
  double e1[2], e2[2], dp[2];
  double d1, d2, x;
  int i;

  e1[0]=sin(azi);
  e1[1]=cos(azi);
  for (i=0; i<2; i++)
    {
      e2[i]=p2[i]-p1[i];
      dp[i]=p1[i]-origin[i];
    }

  x = e1[0]*e2[1] - e1[1]*e2[0];
  if (fabs(x) < 1.0e-15 ) return -1;
  d1 = (dp[0]*e2[1] - dp[1]*e2[0])/x;
  d2 = (dp[0]*e1[1] - dp[1]*e1[0])/x;

  if (!(d1>=0 && d2>=0 && d2<=1.0)) return -1;

  if (ipoint != NULL)
    for (i=0; i<2; i++) 
      ipoint[i]=origin[i] + d1*e1[i];

  return 0;
}

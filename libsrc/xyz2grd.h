/*
 *      FILE: xyz2grd
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: create grid from arbitrary ground points
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: xyz2grd.h,v 1.8 2000/12/11 16:09:00 favey Exp $
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
 *   AUTHORS: Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *               10-1999  --  created
 *            
 *
 */
#ifndef XYZ2GRD
#define XYZ2GRD

#include <stdio.h>
/***************************************************************************/

struct blockchains {
  double bxmin;
  double bxmax;
  double bxinc;
  long int bxn;
  double bymin;
  double bymax;
  double byinc;
  long int byn;
  int *count;
  struct Node ***ch; /* doppeltes Array von Node Pointern */
};

/***************************************************************************/
struct blockchains * blockchain_init (const double xmin, const double xmax,
				      const double ymin, const double ymax,
				      const int xn, const int yn);


/***************************************************************************/
void blockchain_free (struct blockchains *ablock);


/***************************************************************************/
int blockpos (struct blockchains *b, const double xp, const double yp,
	      int *col, int *row);

/***************************************************************************/
int blockchain_insert (struct blockchains *b, double newp[3],
			const int col, const int row);


/***************************************************************************
 *  interpolation for point with coordinates (x,y) as weighted 
 *  average over the circle points (weight = 1/sqr(distance))
 */
double circle_ipol (const double gx, const double gy,
		    struct Node *circle, const int pcount, const int skipg,
		    const int exponent, const double eps);


/***************************************************************************/
int xyz2grd (FILE *fpxyz, FILE *fppts, const char *grdfn,
	     const double xmin, const double xmax,
	     const double ymin, const double ymax,
	     const double xinc, const double yinc,
	     const double radius, const int exponent, const double eps,
	     const int xn, const int yn,
	     const int surfer, const int quadrant,
	     const int minpt, const int dodiff,  const int skipg,
	     const int schlauch, const int verbose);




#endif


/****************************************************************************

Routines:	routines necessary to run the programs coneccted with
		the robot calibration project

Copyright:      (c) 1995 Geodesy and Geodynamics Lab,
                     ETH Zurich

Author:		Etienne Favey & Andreas Schlatter
Tutor:		Hans-Gerd Maas

Address:	Institute of Geodesy and Photogrammetry
		ETH - Hoenggerberg
		CH - 8093 Zurich

Creation Date:	11/1995

Description:		
	
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
****************************************************************************/


/* used in closest pair, source: clop.c */
#define DISTANCE(a,b) sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))

typedef struct		/* used in closest pair, source: closestpair.c */
{ 
	long int pnr;
	double x;
	double y;
} coord_2d;

double	clop();		/* closest pair algorhithm, source: closestpair.c */

int countinbox();	/* find in box a set of points, source: findinbox.c */


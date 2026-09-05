/****************************************************************************
 *
 *  File:	$Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/closestpair.c,v $
 *
 *  Routine:	Find Closest Pair
 *
 *  Copyright:      (c) 1995 Geodesy and Geodynamics Lab,
 *                    ETH Zurich
 *
 *  Author:	Etienne Favey & Andreas Schlatter
 *  Tutor:	H.-G. Maas
 *
 *  Address:	Institute of Geodesy and Photogrammetry
 *		ETH - Hoenggerberg
 *		CH - 8093 Zurich
 *
 *  $Id: closestpair.c,v 1.5 1999/09/08 09:32:31 favey Exp $
 *
 *  Creation Date: 11/1995
 *
 *  Project:	VTB Robot Calibration
 *
 *  Description:This routine attempts to find that pair of points
 *		(out of a coordinate list), whose distance between
 *		each other is minimal, not accepting those smaller
 *		than a lower limit (lowlim). It returns the points
 *		in chigh and clow. The function's return value is
 *		the found smallest distance or 99999 standing for
 *		infinity.
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
 * $Log: closestpair.c,v $
 * Revision 1.5  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.4  1999/06/02 16:13:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1998/05/27 09:28:15  favey
 * ported from Sun version of photogrammetry
 *
 * Revision 1.2  1998/05/27 08:43:03  favey
 * *** empty log message ***
 *
 *	
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "closestpair.h"

int compyqueue(const void *a, const void *b)
{
	if ( ((coord_2d *)a)->y > ((coord_2d *)b)->y ) return (1);
	if ( ((coord_2d *)a)->y < ((coord_2d *)b)->y ) return (-1);
	return (0);
}



double clop(
	coord_2d	*yqueue,	/* yqueue */
	int		n,		/* number of points total */
	double		lowlim,		/* lower limit for delta */
	coord_2d	*chigh,		/* closest high */
	coord_2d	*clow)		/* closest low */

{
	coord_2d	*xtable;	/* xtable, typedef in robix.h */
	int	t=0, c=2, tnb=1, cnb;	/* tail, current, last spot, current in table */
	double		delta, newdelta;
	coord_2d	chi, clo;	/* closest high, closest low */
	coord_2d	current, check;
	int		i, j;

	/* sorting array in y */

	qsort((coord_2d *) yqueue, n, sizeof(coord_2d), compyqueue);

	xtable = (coord_2d *) malloc(n*sizeof(coord_2d));

	chi = yqueue[0];
	clo = yqueue[1];
	delta = DISTANCE(chi,clo);	/* makro DISTANCE defined in robix.h */
	if (delta <= lowlim) delta = 99999;  /* 99999 stands for infinit distance */

	/* initialize x-table */
	xtable[1] = xtable[0] = chi;	/* smallest x first */
	cnb = (chi.x < clo.x) ? (xtable[1] = clo, 1): (xtable[0] = clo, 0);

	while (c < n)
	{ /* transition */
		current = yqueue[c++];

		/* insert the new point into the x-table */
		for (i=tnb; ((i>=0) && (xtable[i].x > current.x)); i--)
			xtable[i+1] = xtable[i];
		xtable[(cnb = i+1)] = current;
		tnb++;

		/* remove points outside the delta-slice from the x-table */
		while (delta <= current.y - yqueue[t].y)
		{
			for (i=0; (xtable[i].pnr != yqueue[t].pnr); i++) 
			{} /* find tail in xtable */

			for (j=i; (j<tnb); j++) 
			{	
				xtable[j] = xtable[j+1];
				if ((j+1) == cnb) cnb--;
			}
			t++; tnb--;
		}

		/* check the successors of the new point in the x-table */
		if ((i=cnb)<tnb)
			do
			{
				check = xtable[++i];
				newdelta = DISTANCE(current, check);
				if ((newdelta < delta) && (newdelta > lowlim))
				{
					delta = newdelta;
					chi = check;
					clo = current;
				}
			} while ((check.x - current.x > delta) && (i<tnb));
		
		/* check the predecessors of the new point in the y-table */
		if ((i=cnb)>0)
			do
			{
				check = xtable[--i];
				newdelta = DISTANCE(current, check);
				if ((newdelta < delta) && (newdelta > lowlim))
				{
					delta = newdelta;
					chi = check;
					clo = current;
				}
			} while ((current.x - check.x > delta) && (i>0));
	}

	free(xtable);
	*chigh = chi; *clow = clo;
	return(delta);
}




/****************************************************************************
 *
 *  Routine:		Find Points in Box
 *
 *  Author:		Etienne Favey & Andreas Schlatter
 *  Tutor:		H.-G. Maas
 *
 *  Address:		Institute of Geodesy and Photogrammetry
 *			ETH - Hoenggerberg
 *			CH - 8093 Zurich
 *
 *  Creation Date:	11/1995
 *
 *  Project:		VTB Robot calibration
 *
 *  Description:	Given a set of point coordinates and one of those as
 *			'hook', this procedure looks for 4 other points in
 *			in a quadratic box of size boxsize. Look at the
 *			if sequence to see where this box is in relation to
 *			hook.
 *	
 ****************************************************************************/

int countinbox(
	coord_2d	*plist,		/* list of points */
	int		inlist,		/* number of points in list */
	int		n,		/* number of points in inbox */
	int		boxsize,	/* boxsize, you fool */
	coord_2d	hook,		/* point defining box position */
	coord_2d	*found,		/* found points will be returned */
	FILE		*fplog)		/* log file for debugging purposes */
{
  int	
    c=0, i=0;

  while ((c < n) && (i < inlist))
    {
      if    ((   plist[i].y           > hook.y)
	     && (plist[i].y-boxsize   < hook.y)
	     && (plist[i].x+boxsize/2 > hook.x)
	     && (plist[i].x-boxsize/2 < hook.x))
	{
	  c++;
	  found[c]=plist[i];
#ifdef WITH_DEBUG
	  if (fplog) fprintf(fplog,"found %d: %ld\n",c,plist[i].pnr);
#endif
	}
      i++;
    }
  return (c == n);
}

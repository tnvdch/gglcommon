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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "robix.h"

int countinbox(
	coord_2d	*plist,		/* list of points */
	int		inlist,		/* number of points in list */
	int		n,		/* number of points in inbox */
	int		boxsize,	/* boxsize, you fool */
	coord_2d	hook,		/* point defining box position */
	coord_2d	*found,		/* found points will be returned */
	FILE		*fplog)		/* log file for debugging purposes */
{
	int	c=0, i=0;

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
			if (fplog) fprintf(fplog,"found %d: %d\n",c,plist[i].pnr);
#endif
		}
		i++;
	}
	return (c == n);
}

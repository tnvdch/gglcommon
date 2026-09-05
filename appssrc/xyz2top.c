/*
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *  AUTHOR: Marc Cocard <cocard@geod.ethz.ch>
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <common.h>

void main ()
{
	double t, d, x[3], n[3];
	while (fscanf(stdin," %lf %lf %lf %lf %lf", 
		      &t,&(x[0]),&(x[1]),&(x[2]),&d)!=EOF) {
		xyz_to_top (
	 47.4058655681/180.0*M_PI,
	  8.6332130105/180.0*M_PI,
		 x, n);
		fprintf (stdout, "%8.2f %10.6f %10.6f %10.6f %10.6f\n",
			 t, n[0], n[1], n[2], d);
	}
}

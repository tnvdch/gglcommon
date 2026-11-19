/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/d1linear.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: Linear one dimensional interpolation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: d1linear.h,v 1.3 1999/09/08 09:32:31 favey Exp favey $
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            28-04-1998  --  created
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/d1linear.h $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.3  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.2  1998/04/28 16:02:11  favey
 * *** empty log message ***
 *
 *
 *
 */

#ifndef COMMON_D1LINEAR
#define COMMON_D1LINEAR

/*
 *  Simple linear interpolation. 
 */
int d1ipol_linear (
	const double x1,
	const double x2,
	const double y1,
	const double y2,
	const double xip,
	double *yres)
;


#endif

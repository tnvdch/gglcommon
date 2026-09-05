/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapcut.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: fill target map from L+T pixel maps
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapcut.h,v 1.4 2000/08/15 16:38:13 favey Exp $
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: mapcut.h,v $
 * Revision 1.4  2000/08/15 16:38:13  favey
 * mapras: bugfix and complete rewrite
 *
 * Revision 1.3  2000/03/06 17:01:42  favey
 * GMT colortable and grid output
 *
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */

/* test if all sheets are available in the specified directory */
int maptest (const char *mappath,const int anz_karten,int *sheets);

/* calculate the size of the variable outmap */
int map_size(
	const double west,const double north,
	const double east,const double south,
	const double spacing,
	int *outmap_length,int *outmap_width);

void rrborder (
	int sheet_nr, const int map_type,
	const double west,const double north,
	const double east,const double south,
	double *rrwest, double *rrnorth,
	double *rrsouth, double *rreast,
	double *sheet_west, double *sheet_north);

/*zusammenfuegen der einzelnen Ausschnitte der jeweiligen Karten */
int mappaste (
	const double west,const double north,
	const double east,const double south,
	const char *mappath, const char *outfname,
	const int map_type, const int do_gmt, 
	const int lpmm, const int verbose,
	const int anz_karten, int *sheets);


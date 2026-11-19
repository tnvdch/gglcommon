/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapnr.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: determine map number from swiss coordinates
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapnr.h,v 1.2 1999/09/08 09:12:18 favey Exp favey $
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
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/mapnr.h $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */
#define BREITE25 (17500.0)
#define HOEHE25 (12000.0)
#define BREITE100 (4*17500.0)
#define HOEHE100 (4*12000.0)

int nrsheets25 (const double west,const double north,
		const double east,const double south,
		int *sheets, const int size);
void edge_25 (int sheet_nr, double *links, double *oben);
int map_25 (const double y, const double x) ;
int nrsheets100 (const double west,const double north,
		 const double east,const double south,
		 int *sheets, const int size);
void edge_100 (int sheet_nr, double *links, double *oben);
int map_100 (const double y, const double x) ;

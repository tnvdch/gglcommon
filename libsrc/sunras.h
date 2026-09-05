/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sunras.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: I/O of SUN raster graphics files
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sunras.h,v 1.4 2001/02/21 14:00:19 cocard Exp $
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
 * $Log: sunras.h,v $
 * Revision 1.4  2001/02/21 14:00:19  cocard
 * *** empty log message ***
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
struct rasterfile {
  unsigned int ras_magic;              /* magic number */
  int     ras_width;              /* width (pixels) of image */
  int     ras_height;             /* height (pixels) of image */
  int     ras_depth;              /* depth (1, 8, or 24 bits) of pixel */
  int     ras_length;             /* length (bytes) of image */
  int     ras_type;               /* type of file; see RT_* below */
  int     ras_maptype;            /* type of colormap; see RMT_* below */
  int     ras_maplength;          /* length (bytes) of following map */
  /* color map follows for ras_maplength bytes, followed by image */
};

#define RT_BYTE_ENCODED 2       /* Run-length compression of bytes */
#define RMT_NONE        0       /* ras_maplength is expected to be 0 */

int read_sunras_header (FILE *fpras, struct rasterfile *header);
int write_sunras_header (FILE *fpras, struct rasterfile *header, const int inv);

int read_colormap (FILE *fpras, unsigned int colmapsize,  unsigned char *colmap);
int write_colormap (FILE *fpras,unsigned int colmapsize,  unsigned char *colmap);
int write_cptfile (FILE *fpcpt, unsigned int colmapsize, unsigned char *colmap);

void write_sunras_info ( struct rasterfile rh);

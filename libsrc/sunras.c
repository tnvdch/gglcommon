/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sunras.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: I/O of SUN raster graphics files
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sunras.c,v 1.4 2001/02/21 14:00:19 cocard Exp $
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
 * $Log: sunras.c,v $
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
/* erstellt im SS99 waehrend des Vertiefungsblockes Flottenmanagement */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errproc.h"
#include "sunras.h"

#define RAS_MAGIC     0x59a66a95
#define INV_RAS_MAGIC 0x956aa659
#define RT_BYTE_ENCODED 2       /* Run-length compression of bytes */

void write_sunras_info ( struct rasterfile rh)
{
	fprintf(stderr, "%d : ras_magic\n", rh.ras_magic);
	fprintf(stderr, "%d : ras_width\n", rh.ras_width);
	fprintf(stderr, "%d : ras_height\n", rh.ras_height);
	fprintf(stderr, "%d : ras_depth\n", rh.ras_depth);
	fprintf(stderr, "%d : ras_length\n", rh.ras_length);
	fprintf(stderr, "%d : ras_type\n", rh.ras_type);
	fprintf(stderr, "%d : ras_maptype\n", rh.ras_maptype);
	fprintf(stderr, "%d : ras_maplength\n", rh.ras_maplength);
}

static void invert_header (struct rasterfile *header)
{
        int i,j;
        unsigned char *tmp=(unsigned char *)header;
        unsigned char b;
        int itype = sizeof(int);

        for (i=0;i<8;i++)
                for (j=0;j<itype/2;j++)
                { 
                        b = tmp[i*itype+j];
                        tmp[i*itype+j] = tmp[i*itype+itype-j-1];
                        tmp[i*itype+itype-j-1] = b;
                }
}


int read_sunras_header (FILE *fpras, struct rasterfile *header)
{
	int inv=0;

	if (fread ((void *)header, sizeof (struct rasterfile), 
		   (size_t)1, fpras) != 1) {
		puterr ("sunras", "read_sunras_header");
		fprintf (stderr, "I/O Error!\n");
		return 2;
	}
	if (header->ras_magic == INV_RAS_MAGIC) {
		invert_header(header);
		inv=1;
	}
	if (header->ras_magic != RAS_MAGIC) {
		puterr ("sunras", "read_sunras_header");
		fprintf (stderr, "invalid magic number in SUN rasterfile!\n");
		return 2;
	}
	if (header->ras_type == RT_BYTE_ENCODED) {
		puterr ("sunras", "read_sunras_header");
		fprintf (stderr, "run-length-encoded sunras not supported!\n");
		return 2;
	}
	return inv;
}

int write_sunras_header (FILE *fpras, struct rasterfile *header, const int inv)
{
	if (inv == 1) {
		invert_header(header);
	}
	if (fwrite ((void *)header, sizeof (struct rasterfile), 
		    (size_t)1, fpras) != 1) {
		puterr ("sunras", "write_sunras_header");
		fprintf (stderr, "I/O Error!\n");
		return 2;
	}
	if (inv == 1) {
		invert_header(header);
	}
	return 0;
}

int read_colormap (FILE *fpras, unsigned int colmapsize,  unsigned char *colmap)
{
	if (fread ((void *)colmap, sizeof (unsigned char), 
		   (size_t)colmapsize, fpras) != colmapsize) {
		puterr ("sunras", "read_colormap");
		fprintf (stderr, "I/O Error!\n");
		return 2;
	}
	return 0;
}

int write_colormap (FILE *fpras, unsigned int colmapsize,  unsigned char *colmap)
{
	if (fwrite ((void *)colmap, sizeof (unsigned char), 
		    (size_t)colmapsize, fpras) != colmapsize)
	{
		puterr ("sunras", "write_colormap");
		fprintf (stderr, "I/O Error!\n");
		return 2;
	}
	return 0;
}

int write_cptfile (FILE *fpcpt, unsigned int colmapsize, unsigned char *colmap)
{
  extern char *program_invocation_short_name;
  int i;
  if (colmapsize!=768) {
    puterr ("sunras","write_cptfile");
    fprintf (stderr, "Unknown colormap type");
    return 2;
  }
  fprintf (fpcpt,"#\tcpt file crated by %s\n", program_invocation_short_name);
  fprintf (fpcpt,"#\tSun Raster Colormap was given\n#\n");
  for (i=1;i<256;i++) {
    fprintf (fpcpt,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	     i-1, *(colmap+i-1), *(colmap+256+i-1), *(colmap+512+i-1),
	     i,   *(colmap+i),   *(colmap+256+i),   *(colmap+512+i)   );
  }
  fprintf (fpcpt,"B\t0\t0\t0\nF\t0\t0\t0\n");
  return 0;
}

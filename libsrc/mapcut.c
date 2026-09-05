/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapcut.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: fill target map from L+T pixel maps
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapcut.c,v 1.4 2000/08/15 16:38:13 favey Exp $
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
 * $Log: mapcut.c,v $
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
/* erstellt im SS99 waehrend des Vertiefungsblockes Flottenmanagement */

/* Width and height of images in pixels */
#define P_WIDTH (700)
#define P_HEIGHT (480)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "errproc.h"
#include "inout.h"
#include "sunras.h"
#include "mapnr.h"
#include "mapcut.h"

int maptest (const char *mappath,
	     const int anz_karten, int *sheets)
     /* test if all sheets are available in the specified directory */
{
  char fname[1000];
  FILE *fp;
  int lv, ispipe;
  for (lv=0;lv<anz_karten;lv++) {
    sprintf (fname,"%s/komb%d.ras",mappath,sheets[lv]);
    fp = openzipped(fname,&ispipe);
    if (fp==NULL) return -1;
    if (ispipe) { pclose (fp); } else { fclose (fp); }
  }
  return 0;
}


int map_size(const double west,const double north,
	     const double east,const double south,
	     const double spacing,
	     int *outmap_height, int *outmap_width)
/* calculate the size of the variable outmap */
{
  *outmap_width  = (int) floor((east-west)/spacing);
  *outmap_height = (int) floor((north-south)/spacing);
  return (*outmap_width) * (*outmap_height);
}



void rrborder (int sheet_nr, const int map_type,
	       const double west,const double north,
	       const double east,const double south,
	       double *rrwest, double *rrnorth,
	       double *rreast, double *rrsouth,
	       double *sheet_west, double *sheet_north)
{

  switch (map_type) {
  case 25:
    edge_25 (sheet_nr, sheet_west, sheet_north);
    break;
  case 100:
    edge_100 (sheet_nr, sheet_west, sheet_north);
    break;
  }
  if ( *sheet_west < west) {
    *rrwest=west;
  } else {
    *rrwest=*sheet_west;
  }
  if ( *sheet_north > north) {
    *rrnorth=north;
  } else {
    *rrnorth=*sheet_north;
  }
  if ( (*sheet_west+(P_WIDTH)*map_type) < east) {
    *rreast=(*sheet_west+(P_WIDTH)*map_type);
  } else {
    *rreast=east;
  }
  if ( *sheet_north-(P_HEIGHT)*map_type < south) {
    *rrsouth=south;
  } else {
    *rrsouth=(*sheet_north-(P_HEIGHT)*map_type);
  }
}


static int overread_chars (FILE *fp, const int count)
{
  int i, n;
  for (i=0;i<count;i++) {
    n=getc(fp);
    if (n== EOF) return -1;
  }
  return 0;
}

int mappaste (const double west,const double north,
	      const double east,const double south,
	      const char *mappath, const char *outfname,
	      const int map_type, const int do_gmt,
	      const int lpmm, const int verbose,
	      const int anz_karten, int *sheets)
	/*zusammenfuegen der einzelnen Ausschnitte der jeweiligen Karten */
{
  FILE *fp,*fpout;
  int lv, zip, err;
  char fname[246];  
  unsigned char colmap[800];
  unsigned char *outmap;
  double N, E, rrwest, rrnorth, rreast, rrsouth;
  int outmap_height, outmap_width, outmap_size;
  double sheet_west, sheet_north;
  long int n;
  int inv=0;
  struct rasterfile rh;  
  long int pic, pir, poc, por, pip, pop, curpos, toread;
  double spacing = (double) map_type / (double) lpmm;

  outmap_size=map_size(west, north, east, south, spacing,
		       &outmap_height, &outmap_width);
  outmap=(unsigned char *) malloc(outmap_size*sizeof(unsigned char));
  if (outmap == NULL) {
    puterr("mapcut","mappaste");
    fprintf (stderr, "Not enough memory available, I need %d kB!\n",
	     outmap_size/1024);   
    return GGLERR_NOMEM;
  }
  for (lv=0;lv<anz_karten;lv++) {
    if (verbose) fprintf(stderr, "processing map %d\n", sheets[lv]);
    /*opening of the mapsheet under consideration if file is zipped or not */
    sprintf (fname, "%s/komb%d.ras", mappath, sheets[lv]);
    fp = openzipped (fname, &zip);
    if (!fp) {
      puterr("mapcut","mappaste");
      fprintf (stderr, "Error opening file %s!\n", fname);
      return EXIT_FAILURE;
    }
    
    /* read the header */
    inv = read_sunras_header (fp, &rh);
    if (inv == 2) return GGLERR_INOUT;
    err = read_colormap (fp, rh.ras_maplength, colmap);
    if (err != 0) return GGLERR_INOUT;
		
    /* read image pixels */
    if (verbose) fprintf(stderr, "reading image pixels\n");
    rrborder(sheets[lv], map_type, west, north, east, south,
	     &rrwest, &rrnorth, &rreast, &rrsouth, &sheet_west, &sheet_north);
    curpos = 0;

    for (N=rrnorth;N>rrsouth;N-=spacing) {
      pir = (int) floor((sheet_north-N)/spacing);
      por = (int) floor((north-N)/spacing);
      for (E=rrwest;E<rreast;E+=spacing) {
	pic = (int) floor((E-sheet_west)/spacing);
	poc = (int) floor((E-west)/spacing);
	pip = pir * P_WIDTH*lpmm + pic;
	pop = por * outmap_width + poc;
	if (pip<0 || pip >= P_WIDTH*lpmm*P_HEIGHT*lpmm || pop<0 || pop >= outmap_size) {
	  puterr("mapcut","mappaste");
	  fprintf (stderr, "internal pixel accounting error!\n");
	  return EXIT_FAILURE;
	}
	toread = pip-curpos;
	if (toread > 0 && overread_chars(fp, toread)!=0) return -1;
	curpos += toread;
	n = getc(fp);
	if (n==EOF) return EXIT_FAILURE;
	curpos++;
	outmap[pop] = (unsigned char) n;
      }
    }

    if (verbose) fprintf(stderr, "map processing finshed\n");
    if (zip) {
      pclose(fp);
    } else {
      fclose(fp);
    }
  } /* for lv */


  if (verbose) fprintf(stderr, "reads the header\n");
  rh.ras_width = outmap_width;
  rh.ras_height = outmap_height;
  rh.ras_length = rh.ras_width*rh.ras_height;
  if (do_gmt == FALSE) {
    /* write sun raster file (MT) original code */
    if (!(fpout = open_text_file (outfname, "w"))) return EXIT_FAILURE;
    write_sunras_header (fpout, &rh, inv);
    write_colormap (fpout, rh.ras_maplength, colmap);
    if (verbose) fprintf(stderr, "%d : map size\n", outmap_size);
    fwrite(outmap, 1, outmap_size, fpout);
    fclose (fpout);
  } else {
    /* write GMT Colortable and grid for grdview... (EF) */
    char *cptfname, *grdfname;
    double easting, northing;
    int i,j;
    cptfname = add_extension (outfname, "cpt");
    if (!(fpout = open_text_file (cptfname, "w"))) return EXIT_FAILURE;
    write_cptfile (fpout, rh.ras_maplength, colmap);
    fclose (fpout);
    grdfname = add_extension (outfname, "xyz");
    if (!(fpout = open_text_file (grdfname, "w"))) return EXIT_FAILURE;
    fprintf (fpout,"# mapras: %9.2f/%9.2f/%9.2f/%9.2f %d\n", 
	     west, east, south, north, outmap_size);
    for (i=0;i<outmap_height;i++)
      for (j=0;j<outmap_width;j++) {
	easting = west+j*1.25;
	northing = north-i*1.25;
	fprintf (fpout, "%9.2f %9.2f %d\n", easting, northing, 
		 *(outmap+(i*outmap_width)+j) );
      }
    fclose (fpout);
  }
  free (outmap);
  if (verbose) fprintf(stderr, "program mapras finished\n");
  return 0;
  
}


/*
 *      FILE: gridop
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: grid operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *            31-08-1997  --  created
 *
 */

#ifndef GGL_GRIDOP
#define GGL_GRIDOP

/* this is for making it compile on Solaris */
#ifndef NAN
#define NAN HUGE_VAL
#endif

/* environment variables for GMT rasters */
#define GMT_RASTERINFO_FILE "RASTERINFO"
#define GMT_GRID_DIR "GMT_GRIDDIR"

/* surfer specific */
#define SURFERMASK (1.70141e38)
#define SURFERMASK_STR " 1.70141e38"

struct Grid {
  char *grdfilename;     /* name of grid file */
  FILE *fp;              /* file pointer to gridfile */
  char *id;              /* identification string */
  long int nx, ny;       /* number of columns and rows */
  double xmin, xmax;     /* x-range */
  double ymin, ymax;     /* y-range */
  double zmin, zmax;     /* y-range */
  double dx, dy;         /* spacing */
  
  double offset, scale;  /* zvalue = offset + data*scale */
  double mask;           /* mask value */

  char datatype;         /* u = unsigned char */
			 /* c = char */
			 /* d = unsigned short int */
			 /* i = short int */
			 /* f = float */
			 /* e = double */      

  unsigned char datasize; /* size of corresponding datatype */
  long int d_offset;      /* offset to first element */
  long int d_row;         /* length of row */
  char *data;             /* pointer to data buffer */
};

/*--------------------------------------------------------------------------------- 
 *  GRID : initialize basic griddata structure without memory allocation for data 
 */
struct Grid * g_init0();

/*--------------------------------------------------------------------------------- 
 *  GRID : initialize memory for griddata depending on datatype  
 */
void * g_getmem4data(struct Grid * g);

/* ---------------------------------------------------------------------------------
 * GRID : set d_offset and d_row according to modus
 *   modus 0 = upside up   --> first element in databuffer = data(xmin, ymin)
 *   modus 1 = upside down --> first element in databuffer = data(xmin, ymax)
 */
void g_setdataparameter(struct Grid *g, int modus);

/*--------------------------------------------------------------------------------- 
 *  GRID : free grid 
 */
void g_free(struct Grid * g);

/*--------------------------------------------------------------------------------- 
 *  GRID : return offset in bytes for grid point (ix, iy)  
 */
long int g_index(struct Grid *g, const long int ix, const long int iy);

/*--------------------------------------------------------------------------------- 
 *  GRID : write grid info to file (for debugging purposes)
 */
void g_printinfo(FILE *fp, struct Grid *g);

/* ---------------------------------------------------------------------------------
 * GRID : scale limits  by factor scalex, scaley resp.
 */
void g_scalelimits(struct Grid *g, double xscale, double yscale);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getval(struct Grid *g, long int offset);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getvalue(struct Grid *g, const long int ix, const long int iy);

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 */
int g_putval(struct Grid *g, long int offset, const double val);

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 */
int g_putvalue(struct Grid *g, const long int ix, const long int iy, const double v);

/* ---------------------------------------------------------------------------------
 * GRID : looks for minimal and maximal values which are stored in g->zmin, g->zmax
 */
void g_minmax(struct Grid *g);

/*---------------------------------------------------------------------------------
 *  GRID : return the 4 corner values and x- and y-values for point 0
 *         ordering of the      y
 *         points 0,1,2,3 :     |
 *                              3 - 2
 *                              | x |
 *                              0 - 1 -> x
 */
int g_getcorner(struct Grid *g, const double x, const double y, 
		double z[4], double * x0, double *y0);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  a vector xyz[3] with x-, y-, and z-value from grid at position (ix,iy) 
*/
int g_getpoint (struct Grid *g, const long int ix, 
                                const long int iy, double xyz[3]);

/* ---------------------------------------------------------------------------------
 * GRID : return the maximal value of the four neighbor
 * only 1 grid point needed
 */
int g_maxof4neighbor(struct Grid *g, const double x, const double y, double *zvalue);

/* ---------------------------------------------------------------------------------
 * GRID : interpolation for point with coordinates (x,y) as weighted average over the
 *        4 corner points (weight = 1/sqr(distance))
 *
 *        g_inter_simple : new version with additional parameter
 *             minpoint = minimal number of corners needed (ranging from 1 - 4)
 */
int g_inter_simple(struct Grid *g, const double x, const double y, int minpoint, double *zvalue);

/*---------------------------------------------------------------------------------
 *  GRID : bilinear interpolation for point with coordinates (x,y)
 */ 
int g_inter_bilinear (struct Grid *g, const double x, const double y, double *zvalue);

/*---------------------------------------------------------------------------------
 *  GMT : returns pointer to start of next substring and length in *n
 *       (may be moved to inout of libcommon)
 */
char * next_substring(char *line, int *n);

/*---------------------------------------------------------------------------------
 *  GMT : decodes a R-option string  (example : -R0/360/-90/90)
 */
int R_option( char * s, double *xmin, double *xmax, double *ymin, double *ymax);


/* ---------------------------------------------------------------------------------
 *  GMT : decodes an I-option string  (example : -I15m)
 */
int I_option( char * s, double *dx, double *dy);

/* ---------------------------------------------------------------------------------
 *  GMT : get info from grdraster.info file
 */
struct Grid * gmt_rasterinfo(int grdnr);

/* ---------------------------------------------------------------------------------
 *  GMT : open rasterfile for reading
 */
int gmt_openrasterfile(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : close rasterfile
 */
void gmt_closerasterfile(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : get memory for data and read whole raster data from file
 */
int gmt_readgrid(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : initialize subgrid and read data from file
 */
struct Grid * gmt_readsubgrid(struct Grid *oldgrid, 
			      double xmin, double xmax, double ymin, double ymax );

/* ---------------------------------------------------------------------------------
 *  SURFER : read surfer grid file ASCII (grid is initialized in routine)
 */
struct Grid *surf_readgrid(const char *grdfname);

/* ---------------------------------------------------------------------------------
 *  SURFER : write surfer grid to file ASCII (format string of 1 number i.e."10.3lf")
 */
int surf_writegrid(const char *grdfname, struct Grid *g, const char *format);

/***************************************************************************/
/*  GMT: write GMT grid to file (uses xyz2grd of GMT)
 *       unit: will be set for all three elements in -D option of xyz2grd
 */
int gmt_writegrid (const char *grdfn, struct Grid *g,
		   const char *unit, const char *comment, const int verbose);

/***************************************************************************/
/*   write grid to ascii file
 *   format entire line with format string. e.g. "%10.3f %10.3f %10.3f\n"
 */
int g_writeasciigrid (const char *grdfn, struct Grid *g, const char *format);

#endif

/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/gridop.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: grid operations (new version sept.98)
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: gridop.c,v 1.22 2002/01/30 13:36:47 cocard Exp $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *            Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            23-11-1997  --  created
 * $Log: gridop.c,v $
 * Revision 1.22  2002/01/30 13:36:47  cocard
 * *** empty log message ***
 *
 * Revision 1.21  2001/07/25 16:20:53  cocard
 * *** empty log message ***
 *
 * Revision 1.20  2001/06/08 12:32:23  cocard
 * *** empty log message ***
 *
 * Revision 1.19  2001/02/21 14:00:19  cocard
 * *** empty log message ***
 *
 * Revision 1.18  2000/08/17 09:53:08  favey
 * gmt_writegrid: allow comment
 *
 * Revision 1.17  2000/03/17 16:06:43  favey
 * new function g_writeasciigrid
 *
 * Revision 1.16  1999/10/11 15:44:33  favey
 * new function gmt_writegrid
 *
 * Revision 1.15  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.14  1999/08/06 16:08:24  cocard
 * *** empty log message ***
 *
 * Revision 1.13  1999/06/10 15:30:31  cocard
 * error in mask checking debugged (mc)
 *
 * Revision 1.12  1999/06/10 11:53:51  cocard
 * *** empty log message ***
 *
 * Revision 1.11  1999/06/02 16:13:19  cocard
 * *** empty log message ***
 *
 * Revision 1.10  1999/02/10 18:03:59  favey
 * changes to make commonlib compile on Solaris
 *
 * Revision 1.9  1999/01/25 16:03:55  favey
 * added put value routines
 *
 * Revision 1.8  1998/12/21 17:29:53  cocard
 * gmt_readsubgrid debugged
 *
 * Revision 1.7  1998/12/21 15:55:22  cocard
 * *** empty log message ***
 *
 * Revision 1.6  1998/12/21 10:06:56  cocard
 * problems with NAN on molnya
 *
 * Revision 1.5  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.4  1998/09/28 15:08:02  cocard
 * new version of gridop (MC)
 *
 * Revision 1.3  1998/08/04 07:00:02  favey
 * *** empty log message ***
 *
 * Revision 1.2  1998/07/21 13:23:52  favey
 * generalized g_interpol_vector
 *
 * Revision 1.1  1998/04/22 12:32:36  favey
 * Initial revision
 *
 * Revision 0.5  1998/01/20 13:16:35  favey
 * removed netcdf stuff from gridop.c
 *
 * Revision 0.4  1998/01/20 10:17:13  favey
 * fixed singularity bug in g_inter_simple()
 *
 * Revision 0.3  1997/12/24 12:10:42  favey
 * followup of V. 0.29
 *
 */

/*
 *  Standard Libraries
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
 *  our special headers
 */

#include "global.h"
#include "errproc.h"
#include "inout.h"
#include "gridop.h"


/* ---------------------------------------------------------------------------------
 * GRID : set d_offset and d_row according to modus
 *   modus 0 = upside up   --> first element in databuffer = data(xmin, ymin)
 *   modus 1 = upside down --> first element in databuffer = data(xmin, ymax)
 */
void g_setdataparameter(struct Grid *g, int modus)
{
  if (modus == 1) 
    {
      g->d_offset = (g->ny - 1)*g->nx;
      g->d_row = -g->nx;	
    }
  else
    {
      g->d_offset = 0;
      g->d_row = g->nx;
    }
  return;
}


/* ---------------------------------------------------------------------------------
 * GRID : evaluate datasize (local function/not in h-file)
 */
int g_datasize(char datatype)
{
  switch (datatype)
    {
    case 'u': return sizeof(unsigned char);
    case 'c': return sizeof(char);
    case 'd': return sizeof(unsigned short int);
    case 'i': return sizeof(short int);
    case 'f': return sizeof(float);
    case 'e': return sizeof(double);
    default : return 0;
   }
}

/* ---------------------------------------------------------------------------------
 * GRID : initialize basic griddata structure without memory allocation for data
 */
struct Grid * g_init0()
{
  struct Grid * g;

  if ((g = (struct Grid *) malloc(sizeof(struct Grid))) == NULL) return NULL;
  memset(g, 0, sizeof(struct Grid));
  return g;
}

/* ---------------------------------------------------------------------------------
 * GRID : initialize memory for griddata depending on datatype
 */
void * g_getmem4data(struct Grid * g)
{
  long int size;

  if (g==NULL || (g->datasize = g_datasize(g->datatype)) == 0) return NULL;
  size = ((long int) g->datasize) * g->nx * g->ny;
  if ((g->data = (char *) malloc(size)) == NULL) 
    {
      fprintf (stderr, "GLL Common Library -- gridop\n");
      fprintf (stderr, "g_getmem4data: error in memory allocation ");
      fprintf (stderr, "(size = %ld bytes)\n",size);
      return NULL;
    }
  memset(g->data, 0, size);   
  return g->data;
}

/* ---------------------------------------------------------------------------------
 * GRID : free grid
 */
void g_free(struct Grid * g)
{
  if (g==NULL) return;
  if (g->data != NULL) free(g->data);
  if (g->id != NULL) free(g->id);
  if (g->grdfilename != NULL) free(g->grdfilename);
  if (g->fp != NULL) fclose(g->fp);
  free(g);
  return;
}


/* ---------------------------------------------------------------------------------
 * GRID : return offset in bytes for grid point (ix, iy) 
 */
long int g_index(struct Grid *g, const long int ix, const long int iy)
{
  return  (g->d_offset + iy*g->d_row + ix)*g->datasize;
}

/*--------------------------------------------------------------------------------- 
 *  GRID : write grid info to file (for debugging purposes)
 */
void g_printinfo(FILE *fp, struct Grid *g)
{
  fprintf(fp,"gridfilename : %s\n",(g->grdfilename==NULL)?"":g->grdfilename);
  fprintf(fp,"file used    : %s\n",(g->fp == NULL)?"no":"yes");
  fprintf(fp,"data read    : %s\n",(g->data == NULL)?"no":"yes");
  fprintf(fp,"idstring     : %s\n",(g->id ==NULL)?"":g->id);
  fprintf(fp,"datasize     : %ld %ld\n",g->nx, g->ny);
  fprintf(fp,"x-range      : %f %f\n", g->xmin, g->xmax);
  fprintf(fp,"y-range      : %f %f\n", g->ymin, g->ymax);
  if (fabs(g->zmin)<1.0e-15 && fabs(g->zmax)<1.0e-15)
    fprintf(fp,"z-range      : unknown\n");
  else
    fprintf(fp,"z-range      : %f %f\n", g->zmin, g->zmax);
  fprintf(fp,"increment    : %f %f\n", g->dx, g->dy);
  fprintf(fp,"offset       : %f\n", g->offset);
  fprintf(fp,"scaling fac. : %f\n",g->scale);
  fprintf(fp,"mask         : %f\n",g->mask);
  fprintf(fp,"datatype     : %c\n",g->datatype);
  fprintf(fp,"datasize     : %d\n",g->datasize);
  return;
}

/* ---------------------------------------------------------------------------------
 * GRID : scale limits  by factor scalex, scaley resp.
 */
void g_scalelimits(struct Grid *g, double xscale, double yscale)
{
  if (g == NULL) return;
  g->xmin *= xscale;  g->xmax *= xscale;  g->dx *= xscale;
  g->ymin *= yscale;  g->ymax *= yscale;  g->dy *= yscale;
  return;
}

/* ---------------------------------------------------------------------------------
 * GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getval(struct Grid *g, long int offset)
{

  char buffer[16];
  char * pdata;
  double data;

  if (g->data == NULL) /* read data from file */
    {
      if (g->fp == NULL || fseek(g->fp, offset, SEEK_SET) != 0) return NAN;
      if (fread(&buffer, g->datasize, 1, g->fp) != 1) return NAN;
      pdata = buffer;
    }
  else pdata = g->data + offset; /* data already in buffer */

  switch (g->datatype)	
    {
    case 'u': data = (double)(* ((unsigned char *) pdata));         break;
    case 'c': data = (double)(* ((char *) pdata));                  break;
    case 'd': data = (double)(* ((unsigned short int *) pdata));    break;
    case 'i': data = (double)(* ((short int *) pdata));             break;
    case 'f': data = (double)(* ((float *) pdata));                 break;
    case 'e': data = (double)(* ((double *) pdata));                break;
    default : return NAN;
    }  

  if (!isnan(g->mask) && 
      fabs((data-g->mask)/((g->mask==0)?1:g->mask))<1.0e-8) return NAN;

   return g->offset + data*g->scale;
}

/* ---------------------------------------------------------------------------------
 * GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getvalue(struct Grid *g, const long int ix, const long int iy)
{
  return g_getval(g, g_index(g, ix, iy));
}

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 * DEBUG XXX test this routine for integrity!
 */
int g_putval(struct Grid *g, long int offset, const double val)
{
	char buffer[16];
	char * pdata;
	double data;

	data = (val - g->offset) / g->scale;

	if (!isnan(g->mask) && isnan(data)) data = g->mask;

	if (g->data == NULL) pdata = buffer;	/* put data to file */
	else pdata = g->data+offset;		/* data already in buffer */

	switch (g->datatype) {
	case 'u': *((unsigned char *) pdata)     = (unsigned char) rint(data); break;
	case 'c': *((char *) pdata)              = (char) rint(data);          break;
	case 'd': *((unsigned short int *) pdata)= (unsigned short int) rint(data); break;
	case 'i': *((short int *) pdata)         = (short int) rint(data);     break;
	case 'f': *((float *) pdata)             = (float) data;               break;
	case 'e': *((double *) pdata)            = (double) data;              break;
	default : return GGLERR_INVAL;
	}  

	if (g->data == NULL) /* put data to file */
	{
		if (g->fp == NULL || fseek(g->fp, offset, SEEK_SET) != 0) 
			return GGLERR_INVAL;
		if (fwrite(&buffer, g->datasize, 1, g->fp) != 1) return GGLERR_INOUT;
	}
	
	return 0;
}

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 * DEBUG XXX test this routine for integrity!
 */
int g_putvalue(struct Grid *g, const long int ix, const long int iy, const double v)
{
  return g_putval(g, g_index(g, ix, iy), v);
}

/* ---------------------------------------------------------------------------------
 * GRID : looks for minimal and maximal values 
 */
void g_minmax(struct Grid *g)
{
  long int i, j;
  double value;

  g->zmin =  HUGE_VAL;
  g->zmax = -HUGE_VAL;
  for (j=0; j<g->ny; j++) 
    for (i=0; i<g->nx; i++) 
      {
	if (isnan(value = g_getval(g, g_index(g,i,j)) )) continue;
	if (value<g->zmin) g->zmin=value;
	if (value>g->zmax) g->zmax=value;
      }
}

/* ---------------------------------------------------------------------------------
 * GRID : return the 4 corner values and x- and y-values for point 0
 *        ordering of the      y
 *        points 0,1,2,3 :     |
 *                             3 - 2
 *                             | x |
 *                             0 - 1 -> x
 * ATTENTION the return value has been changed
 * former  0 = ok 
 * now     the number of corners not NAN are returned
 */
int g_getcorner(struct Grid *g, const double x, const double y, 
		double z[4], double * x0, double *y0)
{
  long int ix, iy;
  int count, i;

  ix = (long int) ((x-g->xmin)/g->dx);
  iy = (long int) ((y-g->ymin)/g->dy);

  for (i=0; i<4; i++) z[i]=NAN;

  if (ix < 0 || iy < 0 ||     /* here we are outside grid */
      ix > (long int) g->nx-2 || 
      iy > (long int) g->ny-2)       return 0;

  if (x0 != NULL) *x0 = g->xmin + ix*g->dx;
  if (y0 != NULL) *y0 = g->ymin + iy*g->dy;
 
  z[0] = g_getval(g, g_index(g, ix,  iy  ));
  z[1] = g_getval(g, g_index(g, ix+1,iy  ));
  z[2] = g_getval(g, g_index(g, ix+1,iy+1));
  z[3] = g_getval(g, g_index(g, ix,  iy+1));

  count =0;
  for (i=0; i<4; i++) 
    if (! isnan(z[i])) count++;

  return count;
}

/* ---------------------------------------------------------------------------------
 * GRID : return  a vector xyz[3] with x-, y-, and z-value from grid at position (ix,iy)
 */
int g_getpoint (struct Grid *g, 
		const long int ix, 
		const long int iy, 
		double xyz[3])
{
  if ( ix>=g->nx || iy>=g->ny) return GGLERR_INVAL;
  xyz[0] = g->xmin + ix*g->dx;
  xyz[1] = g->ymin + iy*g->dy;
  xyz[2] = g_getval(g, g_index(g, ix, iy));
  return 0;
}

/* ---------------------------------------------------------------------------------
 * GRID : return the maximal value of the four neighbor
 * only 1 grid point needed
 */
int g_maxof4neighbor(struct Grid *g, const double x, const double y, double *zvalue)
{
  double z[4];
  double x0, y0;
  int i;

  *zvalue=NAN;
  g_getcorner(g, x, y, z, &x0, &y0);
  for (i=0; i<4; i++)
    if (!isnan(z[i]) && (isnan(*zvalue) || z[i]> *zvalue))
      *zvalue = z[i];
  return (isnan(*zvalue))?-1:0;
}

/* ---------------------------------------------------------------------------------
 * GRID : interpolation for point with coordinates (x,y) as weighted average over the
 *        4 corner points (weight = 1/sqr(distance))
 *        new version with additional parameter
 *                minpoint = minimal number of corners not NAN
 */
int g_inter_simple (struct Grid *g, const double x, const double y, int minpoint, double *zvalue)
{
  const double eps=1.0e-10;
  double x0, y0;
  double dix, diy;
  double a1, b1, a2, b2;
  double z[4], weight[4], sumweight;
  int i;

  if (g_getcorner(g, x, y, z, &x0, &y0) < minpoint) return GGLERR_INVAL;

  dix = (x-x0)/g->dx;
  diy = (y-y0)/g->dy;
  a1 = dix*dix;  a2=(1-dix)*(1-dix);
  b1 = diy*diy;  b2=(1-diy)*(1-diy);
    
  weight[0] = 1.0/(sqrt(a1 + b1) + eps);
  weight[1] = 1.0/(sqrt(a2 + b1) + eps);
  weight[2] = 1.0/(sqrt(a2 + b2) + eps);
  weight[3] = 1.0/(sqrt(a1 + b2) + eps);

  sumweight=0.0;
  for (i=0; i<4; i++) 
    if (! isnan(z[i])) sumweight += weight[i];

  *zvalue=0;
  for (i=0; i<4; i++) 
    if (! isnan(z[i]))
      *zvalue += weight[i]/sumweight * z[i];
  
  return 0;
}

/* ---------------------------------------------------------------------------------
 * GRID : bilinear interpolation for point with coordinates (x,y)
 */
int g_inter_bilinear (struct Grid *g, const double x, const double y, double *zvalue)
{
  double x0, y0;
  double z[4];
  
  if (g_getcorner(g, x, y, z, &x0, &y0) != 0) return GGLERR_INVAL;

  *zvalue =  z[0]
    + ((z[1]-z[0])/g->dx) * (x-x0) 
    + ((z[3]-z[0])/g->dy) * (y-y0)
    + ((z[2]-z[1]-z[3]+z[0])/g->dx/g->dy) * ((x-x0)*(y-y0));

  return 0;
}

/* ---------------------------------------------------------------------------------
 * GMT : returns pointer to start of next substring and length in *n
 *      (may be moved to inout of libcommon)
 */
char * next_substring(char *line, int *n)
{
  int i=0;
  char *p;
  char idchar;

  /* remove trailing blanks -> position p to first non blank char */
  while (isblank(line[i])) i++; 
  p=line+i;

  /* evaluate length of substring */ 
  *n=0;
  if (isspace(p[0])) return NULL;
  if (p[0] == '"' || p[0]=='\'')   /* quotation mark -> keep string together ... */
    {
      idchar = p[0];
      do {
	(*n)++;
	if (p[*n]=='\0') return NULL;
      } while (p[*n] != idchar);  /* ... until closing quotation mark encountered */
      (*n)++;
    }
  else  /* no quotation mark -> isspace is delimiter */
    while (!isspace(p[*n])) (*n)++;  
  return p;
}

/* ---------------------------------------------------------------------------------
 * GMT : decodes a R-option string  (example : -R0/360/-90/90 or 0/360/-90/90)
 */
int R_option( char * s, double *xmin, double *xmax, double *ymin, double *ymax)
{
  char *text;
  unsigned int i;

 
  if ((text = strdup((s[1]=='R')?s+2:s)) == NULL) return -1;
  
  
  for (i=0; i<strlen(text); i++) if (text[i]=='/') text[i]=' ';
  if (sscanf(text,"%lf %lf %lf %lf", xmin, xmax, ymin, ymax) != 4) return -1;
  if (*xmin>*xmax || *ymin>*ymax) return -1;
  free(text);

  return 0;
}

/* ---------------------------------------------------------------------------------
 * GMT : decodes an I-option string  (example : -I15m or 15m)
 */
int I_option( char * s, double *dx, double *dy)
{
  char *p1, *p2;

  p1=(s[1]=='I')?s+2:s;

  p2=p1; *dx = strtod(p2, &p1); if (p1==p2) return -1;
  if (p1[0] == 'm') {(*dx)/=60;   p1++;}  /* minutes converted to degrees */
  if (p1[0] == 's') {(*dx)/=3600; p1++;}  /* seconds converted to degrees */

  if (p1[0] == '\0') *dy = *dx;
  else
    {
      p2=p1+1; *dy = strtod(p2, &p1); if (p1==p2) return -1;    
      if (p1[0] == 'm') {(*dy)/=60;   p1++;}  /* minutes converted to degrees */
      if (p1[0] == 's') {(*dy)/=3600; p1++;}  /* seconds converted to degrees */
    }
  return 0;
}

/* ---------------------------------------------------------------------------------
 * GMT : get and check env. variable (local function/not in h-file)
 */
char * gmt_getenv(char * name)
{
  char *variable;
  if ((variable = getenv(name)) == NULL)
    {
      fprintf (stderr, "GLL Common Library -- gridop\n");
      fprintf (stderr, "environement variable %s not set\n", name);
      return NULL;      
    }
  return variable;
}


/* ---------------------------------------------------------------------------------
 * GMT : add path to filename (local function/not in h-file)
 */
int gmt_addpath(char *longfilename, char *shortfilename)
{
 char *path;
 int i;

  memset(longfilename, 0, MAXLINELENGTH);
  if (shortfilename[0] != '/')
    {
      if ((path=gmt_getenv(GMT_GRID_DIR)) == NULL) return -1;
      strcpy(longfilename, path);
      i=strlen(longfilename)-1;
      if (longfilename[i] != '/') longfilename[i+1]='/';
    }
  strcat(longfilename, shortfilename);
  return 0;
}


/* ---------------------------------------------------------------------------------
 * GMT : get info from grdraster.info file
 */
struct Grid * gmt_rasterinfo(int grdnr)
{
  FILE *fp;
  char longfilename[MAXLINELENGTH];
  char *shortfilename;
  struct Grid * g;
  char line[MAXLINELENGTH];
  char substr[MAXLINELENGTH];
  char *p;
  int count;
  char registration;

  if ((shortfilename=gmt_getenv(GMT_RASTERINFO_FILE)) == NULL) return NULL;
  if (gmt_addpath(longfilename, shortfilename) != 0) return NULL;

  if ((fp = open_text_file(longfilename, "r")) == NULL) return NULL;
  while (!feof(fp))
    {      
      if (read_line_static(fp, line, MAXLINELENGTH) == 0) return NULL;

      /* field  1: grid number used for identification */
      if (line[0] != '#' && strtol(line, &p, 10) == grdnr)
	{
	  if ((g=g_init0()) == NULL) return NULL;
	  if ((p=next_substring(p,&count)) == NULL) return NULL;

	  /* field  2: grid description */
	  p[count-1]='\0';
	  if ((g->id =strdup(p+1)) == NULL) return NULL;
	  p+=count;
	  
	  /* field 3: zunit skipped */
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
     	  p+=count;
	  
          /* field 4: R option (boundaries) */
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  memset(substr, 0, MAXLINELENGTH );
	  strncpy(substr, p, count);
	  if (R_option(substr,&(g->xmin), &(g->xmax),
		              &(g->ymin), &(g->ymax)) != 0) return NULL;
	  p+=count;

	  
	  /* field 5: I option (interval) */
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  memset(substr, 0, MAXLINELENGTH );
	  strncpy(substr, p, count);
	  if (I_option(substr, &(g->dx), &(g->dy)) != 0) return NULL;
	  p+=count;

          /* field 6: Pixel or grid registration */
          if ((p=next_substring(p,&count)) == NULL) return NULL;
	  registration = p[0];
	  if (!(registration == 'P' || registration == 'G')) return NULL;
	  p+=count;

	  /* field 7: data type */
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  g->datatype = p[0];
	  g->datasize = g_datasize(g->datatype);
	  p+=count;

	  /* field 8: scale */
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  g->scale = strtod(p, NULL);
	  p+=count;

	  /* field 9: offset */ 
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  g->offset = strtod(p, NULL);
	  p+=count;

	  /* field 10: mask value */ 
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  if (strncmp(p,"none",4) == 0) g->mask=NAN;
	  else 
	    {
	      g->mask = strtod(p, NULL);
	      if (fabs(g->mask)<1.0e-14) 
		memset(&(g->mask),0,sizeof(g->mask));
	    }
	  p+=count;

	  /* field 11: filename  */ 
	  if ((p=next_substring(p,&count)) == NULL) return NULL;
	  p[count]='\0';
	  if ((g->grdfilename = strdup(p)) == NULL) return NULL;

	  fclose(fp);

	  /* we change pixel registration to grid registration */  
	  if (registration == 'P')
	    {
	      g->xmin += g->dx/2.0;
	      g->xmax -= g->dx/2.0;
	      g->ymin += g->dy/2.0;
	      g->ymax -= g->dy/2.0;
	    }

	  g->nx = (long int) rint((g->xmax - g->xmin)/g->dx + 1);
	  g->ny = (long int) rint((g->ymax - g->ymin)/g->dy + 1);
	  g_setdataparameter(g, 1);
	  return g;
	}
    }
  return NULL;
}

/* ---------------------------------------------------------------------------------
 * GMT : open rasterfile for reading
 */
int gmt_openrasterfile(struct Grid *g)
{
  char filename[MAXLINELENGTH];

  if (gmt_addpath(filename, g->grdfilename) != 0) return -1;
  if ((g->fp = open_text_file(filename,"r")) == NULL) return -1;  
  return 0;
}

/* ---------------------------------------------------------------------------------
 * GMT : close rasterfile
 */
void gmt_closerasterfile(struct Grid *g)
{
  if (g->fp != NULL) fclose(g->fp);
  g->fp = NULL;
}

/* ---------------------------------------------------------------------------------
 * GMT : get memory for data and read whole raster data from file
 */
int gmt_readgrid(struct Grid *g)
{
  unsigned long int n;

  if (g_getmem4data(g) == NULL) return -1;
  if (g->fp == NULL && gmt_openrasterfile(g) !=0) return -1;
  n = g->nx*g->ny;
  
  if (fread(g->data, g->datasize, n, g->fp) != n) 
    {
      fprintf (stderr, "GLL Common Library -- gridop\n");
      fprintf (stderr, "g_readgrid: error in reading data\n");
      return -1;
    }
  gmt_closerasterfile(g);
  return 0;
}

/* ---------------------------------------------------------------------------------
 * GMT : initialize subgrid and read data from file
 */
struct Grid * gmt_readsubgrid(struct Grid *oldgrid, 
			      double xmin, double xmax, double ymin, double ymax )

{
  FILE *fp;
  struct Grid * g;
  char *data;
  long int i;
  long int position;

  /* check if subgrid limits are within gridlimits */
  if (! (xmin >= oldgrid->xmin && xmax <= oldgrid->xmax &&
	 ymin >= oldgrid->ymin && ymax <= oldgrid->ymax))
      {
	fprintf (stderr, "GLL Common Library -- gridop\n");
	fprintf (stderr, "g_readsubgrid: invalid grid limits\n");
	fprintf (stderr, "x: %f <--(%f %f)--> %f\n",
		 oldgrid->xmin, xmin, xmax, oldgrid->xmax);
	fprintf (stderr, "y: %f <--(%f %f)--> %f\n",
		 oldgrid->ymin, ymin, ymax, oldgrid->ymax);
	return NULL;
      }

  /* initialize new grid */
  if ((g = g_init0()) == NULL) return NULL;
  memcpy(g, oldgrid, sizeof(struct Grid));
  g->grdfilename = NULL;
  g->id = strdup(oldgrid->id);
  g->fp = NULL;

  g->xmin = xmin;   g->xmax = xmax; 
  g->ymin = ymin;   g->ymax = ymax;   
  g->nx = (long int) rint((xmax - xmin)/g->dx) + 1;
  g->ny = (long int) rint((ymax - ymin)/g->dy) + 1;

  /* same modus as original grid */
  g_setdataparameter(g,(oldgrid->d_offset==0)?0:1);
 
  if (g_getmem4data(g) == NULL) return NULL;

  /* read data from file */
  data = g->data + g->d_offset*g->datasize;
  fp=oldgrid->fp;

  position = g_index(oldgrid, 
		     (long int) rint((xmin - oldgrid->xmin)/g->dx), 
		     (long int) rint((ymin - oldgrid->ymin)/g->dy));

  for (i=0;i<g->ny;i++) 
    {
      if (oldgrid->data == NULL)  /* read data from file */
	{
	  if (fseek(fp, position, SEEK_SET) ||
	      fread(data, g->datasize, g->nx, fp) != (unsigned long int) g->nx)
	    {
	      fprintf (stderr, "GLL Common Library -- gridop\n");
	      fprintf (stderr, "g_readsubgrid: error in reading data\n");
	      return NULL;
	    }
	}
      else                        /* copy data from oldgrid */	
	memcpy(data, oldgrid->data + position, g->datasize*g->nx);
      
      data += g->datasize*g->d_row; 
      position += oldgrid->datasize*oldgrid->d_row;  
    }
  
  return g;
}

/*-------------------------------------------------------------------------------
 * SURFER : read surfer grid file ASCII (grid is initialized in this routine)
 */
struct Grid *surf_readgrid(const char *grdfname)
{
  FILE *fp;

  long int i;
  char line[MAXLINELENGTH];
  struct Grid *g;
  char * data;
  
  if ((fp = open_text_file(grdfname, "r")) == NULL) return NULL;
  if (read_line_static(fp, line, MAXLINELENGTH) == 0 || strncmp(line,"DSAA",4) != 0) 
    {
      fprintf (stderr, "GLL Common Library -- gridop\n");
      fprintf (stderr, "g_read: invalid grid format in file: %s", grdfname);
      fprintf (stderr, " (format DSAA expected)\n");
      return NULL;
    }
  /* at this point we assume to have a correct grid format */
  if ((g=g_init0()) == NULL) return NULL;
  g->grdfilename = strdup(grdfname);
  g->id = strdup("DSAA\0");
  g->datatype = 'e';
  g->datasize = g_datasize(g->datatype);
  g->mask = SURFERMASK;
  g->scale = 1;
  g->offset = 0;

  if (fscanf (fp, " %ld %ld", &(g->nx), &(g->ny)) != 2) return NULL;
  if (g_getmem4data(g) == NULL) return NULL;

  if (fscanf (fp, " %lf %lf", &(g->xmin), &(g->xmax)) != 2) return NULL;
  if (fscanf (fp, " %lf %lf", &(g->ymin), &(g->ymax)) != 2) return NULL;
  if (fscanf (fp, " %lf %lf", &(g->zmin), &(g->zmax)) != 2) return NULL;

  g->dx = (g->xmax-g->xmin)/(g->nx-1);
  g->dy = (g->ymax-g->ymin)/(g->ny-1);
  g_setdataparameter(g, 0);

  data=g->data;
  for (i=0; i < g->ny*g->nx; i++)
    {	
      if (fscanf (fp, " %lf", (double*) data) != 1) return NULL;
      data += g->datasize;
    }
  fclose (fp);
  fp=NULL;
  return g;
}


/* ---------------------------------------------------------------------------------
 * SURFER : write surfer grid to file ASCII (format string of 1 number i.e."10.3lf")
 */
int surf_writegrid(const char *grdfname, struct Grid *g, const char *format)
{
  FILE *fp;
  long int i,j;
  double value;

  if (g == NULL) return -1;
  if ((fp = open_text_file(grdfname, "w")) == NULL) return -1;

  fprintf (fp, "DSAA\n");
  fprintf (fp, "%ld %ld\n", g->nx, g->ny);
  fprintf (fp, "%f %f\n", g->xmin, g->xmax);
  fprintf (fp, "%f %f\n", g->ymin, g->ymax);

  fprintf (fp, format,g->zmin);
  fprintf (fp, format,g->zmax);
  fprintf (fp, "\n");

  for (j=0; j<g->ny; j++) 
    {
      for (i=0; i<g->nx; i++) 
	{
	  value = g_getvalue(g,i,j);
	  if (isnan(value)) 
	    fprintf(fp, SURFERMASK_STR);
	  else
	    fprintf(fp, format, value);
	}
      fprintf (fp, "\n");
    }
  fclose (fp);
  return 0;
}

/***************************************************************************/
/*  GMT: write GMT grid to file (uses xyz2grd of GMT)
 *       unit: will be set for all three elements in -D option of xyz2grd
 */
int gmt_writegrid (const char *grdfn, struct Grid *g,
		   const char *unit, const char *comment, const int verbose)
{
  FILE *fpgmt;
  char command[256];
  long int i, j;
  double zvalue;
  extern char *program_invocation_short_name;
  sprintf (command,
	   "xyz2grd -G%s -R%f/%f/%f/%f -I%f/%f -D%s/%s/%s/%f/%f/=/%s -Zd %s",
	   grdfn, g->xmin, g->xmax, g->ymin, g->ymax, g->dx, g->dy,
	   unit, unit, unit, g->scale, g->offset, 
	   comment, (verbose?"-V":""));
  fpgmt = popen (command, "w");
  if (!fpgmt) {
    fprintf (stderr, "%s: unable to open pipe to xyz2grd\n",
	     program_invocation_short_name);
    return GGLERR_INOUT;
  }
  for (j=g->ny-1; j>=0; j--) {
    for (i=0; i<g->nx; i++) {
      zvalue = g_getvalue (g, i, j);
      fwrite (&zvalue, 1, sizeof(double), fpgmt);
    }
  }
  return pclose (fpgmt);
}

/***************************************************************************/
/*   write grid to ascii file, only values != NaN are written
 *   format entire line with format string. e.g. "%10.3f %10.3f %10.3f\n"
 */
int g_writeasciigrid (const char *grdfn, struct Grid *g, const char *format)
{
  FILE *fpasc;
  char command[256];
  long int i, j;
  int err;
  double xvec[3];
  extern char *program_invocation_short_name;

  sprintf (command,
	   "# -G%s -R%f/%f/%f/%f -I%f/%f ; scale %f, offset %f",
	   grdfn, g->xmin, g->xmax, g->ymin, g->ymax, g->dx, g->dy,
	   g->scale, g->offset);
  fpasc = open_text_file (grdfn, "w");
  if (!fpasc) return GGLERR_INOUT;
  fprintf (fpasc, "# grid created by %s\n%s\n", 
	   program_invocation_short_name, command);

  for (j=g->ny-1; j>=0; j--) {
    for (i=0; i<g->nx; i++) {
      err = g_getpoint (g, i, j, xvec);
      if (err != FALSE) return err;
      if (!(isnan(xvec[2]))) 
	fprintf (fpasc, format, xvec[0], xvec[1], xvec[2]);
    }
  }
  fclose (fpasc);
  return 0;
}


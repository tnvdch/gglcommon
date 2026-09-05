/*
 *      FILE: xyz2grd
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: create grid from arbitrary ground points
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: xyz2grd.c,v 1.16 2000/12/11 16:09:00 favey Exp $
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
 *   AUTHORS: Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *               10-1999  --  created
 *            
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "global.h"
#include "errproc.h"
#include "inout.h"
#include "chain.h"
#include "gridop.h"
#include "mat3.h"
#include "xyz2grd.h"

/***************************************************************************/
struct blockchains * blockchain_init (const double xmin, const double xmax,
				      const double ymin, const double ymax,
				      const int xn, const int yn)
{
  struct blockchains *ablock;
  long int i,j;

  ablock = (struct blockchains *) malloc (sizeof(struct blockchains));
  if (ablock==NULL) return NULL;
  ablock->bxinc = (xmax - xmin) / xn;   ablock->byinc = (ymax - ymin) / yn;
  ablock->bxn = xn + 2;                 ablock->byn = yn + 2;
  ablock->bxmin = xmin - ablock->bxinc; ablock->bymin = ymin - ablock->byinc;
  ablock->bxmax = xmax + ablock->bxinc; ablock->bymax = ymax + ablock->byinc;
  ablock->count = (int *)malloc(ablock->bxn*ablock->byn*sizeof(int));
  if (ablock->count == NULL) return NULL;
  ablock->ch = (struct Node ***) malloc (ablock->bxn*sizeof(struct Node **));
  if (ablock->ch == NULL) return NULL;
  for (i=0; i<ablock->bxn; i++) {
    ablock->ch[i]=(struct Node **)malloc(ablock->byn*sizeof(struct Node *));
    if (ablock->ch[i]==NULL) return NULL;
    for (j=0; j<ablock->byn; j++) {
      ablock->ch[i][j] = chain_init ( &(ablock->count[i*ablock->byn+j]) );
      if (ablock->ch[i][j] == NULL) return NULL;
    }
  }
  return ablock;
}

/***************************************************************************/
void blockchain_free (struct blockchains *ablock)
{
  int i,j;

  if (ablock == NULL) return;
  if (ablock->ch != NULL) {
    for (i=0;i<ablock->bxn;i++)
      if (ablock->ch[i] != NULL) {
	for (j=0;j<ablock->byn;j++)
	  if (ablock->count != NULL)
	    chain_free (ablock->ch[i][j],TRUE,ablock->count[i*ablock->byn+j]);
	free (ablock->ch[i]);
      }
    free (ablock->ch);
  }
  if (ablock->count != NULL) free (ablock->count);
  free (ablock);
}


/***************************************************************************/
static void blockchain_info (struct blockchains *ablock, const int verbose)
{
  int i,j;

  fprintf (stderr, "# initialized blockchain info:\n");
  fprintf (stderr, "block-R: %12.3f %12.3f %12.3f %12.3f\n",
	   ablock->bxmin, ablock->bxmax, ablock->bymin, ablock->bymax);
  fprintf (stderr, "block-I: %10.3f %10.3f\n", ablock->bxinc, ablock->byinc);
  fprintf (stderr, "block-b: %ld %ld\n", ablock->bxn, ablock->byn);
  if (verbose) {
    for (j=ablock->byn-1;j>=0;j--) {
      for (i=0;i<ablock->bxn;i++) {
	fprintf (stderr, "%5d ", ablock->count[i*ablock->byn+j]);
      }
      fprintf (stderr, "\n");
    }
  }
}

/***************************************************************************/
int blockpos (struct blockchains *b, const double xp, const double yp,
	      int *col, int *row)
{
  *col = (int) floor ( (xp - b->bxmin) / b->bxinc );
  *row = (int) floor ( (yp - b->bymin) / b->byinc );
  if ( (*col < 0 || *col > b->bxn - 1) ||
       (*row < 0 || *row > b->byn - 1) ) {
    *col = *row = -1;
    return -1;
  }
  return 0;
}


/***************************************************************************/
int blockchain_insert (struct blockchains *b, double newp[3],
			const int col, const int row)
{
  double *newdata;
  struct Node *newnode; 
  newdata = (double *) malloc (3*sizeof(double));
  if (newdata==NULL) return GGLERR_NOMEM;
  v3_copy (newp, newdata);
  newnode = chain_newelement (newdata);
  if (newnode==NULL) { free(newdata); return GGLERR_NOMEM; }
  chain_insert (b->ch[col][row], newnode, &(b->count[col*b->byn+row]));
  return 0;
}


/***************************************************************************/
static struct Grid * grdinit ( const double xmin, const double xmax,
			       const double ymin, const double ymax,
			       const double xinc, const double yinc)
{
  struct Grid *g;
  void *d;

  g = g_init0();
  if (g==NULL) return NULL;
  g->nx = ( (long int) floor ((xmax - xmin) / xinc) ) + 1;
  g->ny = ( (long int) floor ((ymax - ymin) / yinc) ) + 1;
  g->xmin = xmin;     g->xmax = xmin + (g->nx - 1) * xinc;
  g->ymin = ymin;     g->ymax = ymin + (g->ny - 1) * yinc;
  g->dx = xinc; g->dy = yinc;
  g->offset = 0.0; g->scale = 1.0;
  g->mask = 0.0;
  g->datatype = 'e';
  d = g_getmem4data (g);
  if (d == NULL) { g_free(g); g=NULL; }
  g_setdataparameter (g, 0);
  return g;
}

/***************************************************************************/
#define DIST(x1,y1,x2,y2) sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
static int chainup (const double gx, const double gy, const double r,
		     struct Node *circle, int *pcount,
		     struct Node *block)
{
  struct Node *tmp, *newnode;
  tmp = block;
  while ( (tmp = tmp->next)!= NULL ) {
    if (DIST(gx,gy,((double *)(tmp->data))[0],((double *)(tmp->data))[1])<r) {
      newnode = chain_newelement (tmp->data);
      if (newnode == NULL) return GGLERR_NOMEM;
      chain_insert (circle, newnode, pcount);
    }
  }
  return 0;
}
#undef DIST

/***************************************************************************/
static int check_quadrant (const double gx, const double gy,
		    struct Node *circle, const int pcount,
		    const int quad, const int minpt)
{
  struct Node *tmp;
  double pt[3];
  int tl = 0, bl = 0, tr = 0, br = 0;

  if (quad != FALSE) { /* quadrant method */
    tmp = circle;
    /*   fprintf (stdout, "# %12.3f %12.3f\n", gx, gy); */
    while ((tmp = tmp->next) != NULL) {
      v3_copy (((double *)(tmp->data)), pt);
      if (gx - pt[0] > 0.0) { /* right */
	if (gy - pt[1] > 0.0) { /* top */
	  tr++;
	} else { /* bottom */
	  br++;
	}
      } else { /* left */
	if (gy - pt[1] > 0.0) { /* top */
	  tl++;
	} else { /* bottom */
	  bl++;
	}
      }
      /*     fprintf (stdout, "%12.3f %12.3f\n", pt[0], pt[1]); */
    }
    if (tr && br && tl && bl) return TRUE;
  } else { /* min number of points method */
    if (pcount >= minpt) return TRUE;
  }
  return FALSE;
}


/***************************************************************************
 *  interpolation for point with coordinates (x,y) as weighted 
 *  average over the circle points (weight = 1/sqr(distance))
 */
double circle_ipol (const double gx, const double gy,
		    struct Node *circle, const int pcount, const int skipg,
		    const int exponent, const double eps)
{
  struct Node *tmp;
  double dix, diy;
  double dist, zvalue;
  double *weight, *z, sumweight;
  double eps2, dist2;
  int i, j;

  z = (double *) malloc (pcount*sizeof(double));
  weight = (double *) malloc (pcount*sizeof(double));
  tmp = circle; i=0;
  eps2 = eps;
  for (j=1;j<exponent;j++) eps2 = eps2*eps;
  while ( (tmp = tmp->next) != NULL) {
    dix = ((double *)(tmp->data))[0]-gx;
    diy = ((double *)(tmp->data))[1]-gy;
    dist2 = dist = sqrt(dix*dix + diy*diy);
    for (j=1;j<exponent;j++) dist2 = dist2*dist;

    z[i] = ((double *)(tmp->data))[2];
    weight[i] = eps2/(dist2 + eps2);
    if (skipg != FALSE && dist < 0.1) weight[i] = 0.0;
    i++;
  }
  if (i != pcount) { 
    fprintf (stderr, "circle_ipol: internal bookkeeping error!");
    exit (0);
  }
  sumweight=0.0;
  for (i=0; i<pcount; i++) sumweight += weight[i];

  zvalue=0;
  for (i=0; i<pcount; i++) zvalue += weight[i]/sumweight * z[i];
  free (z); free (weight);
  return zvalue;
}


/***************************************************************************/
static long int blockchain_fillgrid (struct blockchains *b,
				     struct Grid *g,
				     const double r,   /* radius  */
				     const int exponent,
				     const double eps, /* epsilon */
				     const int quad,
				     const int minpt)
{
  long int grow, gcol;
  double gx, gy;
  int brow, bcol;
  int pcount;
  long int nancount=0;
  struct Node *gchain;
  int top, bottom, left, right;

  for (gcol=0; gcol<g->nx; gcol++) {
    for (grow=0; grow<g->ny; grow++) {
      gx = g->xmin+gcol*g->dx;
      gy = g->ymin+grow*g->dy;
      pcount = 0;
      gchain = chain_init (&pcount);
      if (gchain==NULL) return GGLERR_NOMEM;
      blockpos (b, gx, gy, &bcol, &brow);
      if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow]) != 0) 
	return GGLERR_NOMEM;
      top    = (gy >  (b->bymin + (brow+1)*b->byinc) - r);
      bottom = (gy <= (b->bymin + (brow  )*b->byinc) + r);
      right  = (gx >  (b->bxmin + (bcol+1)*b->bxinc) - r);
      left   = (gx <= (b->bxmin + (bcol  )*b->bxinc) + r);
      if (top)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow]) != 0)
	  return GGLERR_NOMEM;
      if (left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow]) != 0)
	  return GGLERR_NOMEM;
      if (top && right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (top && left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom && right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom && left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (check_quadrant (gx, gy, gchain, pcount, quad, minpt)) {
	g_putvalue (g, gcol, grow, circle_ipol (gx, gy, gchain, pcount, FALSE, exponent, eps));
      } else {
	g_putvalue (g, gcol, grow, g->mask);
	nancount++;
      }
      chain_free (gchain, FALSE, pcount);
    }
  }

  return -nancount;
}

/***************************************************************************/
static long int blockchain_procpts (const double xmin,
				    const double xmax,
				    const double ymin,
				    const double ymax,
				    struct blockchains *b,
				    FILE *fppts,
				    const char *grdfn,
				    const double r,   /* radius  */
				    const int exponent,
				    const double eps, /* epsilon */
				    const int quad,
				    const int minpt,
				    const int dodiff,
				    const int skipg,
				    long int *totpt)
{
  FILE *fpgrd;
  double gx, gy, gh;
  int brow, bcol;
  int pcount, readcount;
  long int nancount=0;
  long int totcount=0;
  struct Node *gchain;
  int top, bottom, left, right;
  char entireline[1024];

  fpgrd = open_text_file (grdfn, "w");
  if (fpgrd == NULL) return GGLERR_INOUT;
  fprintf (fpgrd, "# gx  gy  result");
  if (dodiff != FALSE) fprintf (fpgrd, "  gh  gh-result");
  fprintf (fpgrd, "\n");

  /* read input file and fill blockchains */
  while (read_line_static (fppts, entireline, 1024) != EOF) {
    if (entireline[0] == '#' || entireline[0] == '$' ||
	entireline[0] == '\n'|| entireline[0] == '\r')
      continue;
    if (dodiff != FALSE) {
      if (sscanf (entireline, " %lf %lf %lf%n", &gx, &gy, &gh, &readcount) < 3)
	return GGLERR_INOUT;
    } else {
      if (sscanf (entireline, " %lf %lf%n", &gx, &gy, &readcount) < 2)
	return GGLERR_INOUT;
    }
    if ( (xmin <= gx) && (gx <= xmax) &&
	 (ymin <= gy) && (gy <= ymax)   ) {
      /* process this point */
      blockpos(b, gx, gy, &bcol, &brow);
      totcount++;
      pcount = 0;
      gchain = chain_init (&pcount);
      if (gchain==NULL) return GGLERR_NOMEM;
      if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow]) != 0) 
	return GGLERR_NOMEM;
      top    = (gy >  (b->bymin + (brow+1)*b->byinc) - r);
      bottom = (gy <= (b->bymin + (brow  )*b->byinc) + r);
      right  = (gx >  (b->bxmin + (bcol+1)*b->bxinc) - r);
      left   = (gx <= (b->bxmin + (bcol  )*b->bxinc) + r);
      if (top)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow]) != 0)
	  return GGLERR_NOMEM;
      if (left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow]) != 0)
	  return GGLERR_NOMEM;
      if (top && right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (top && left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow+1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom && right)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol+1][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (bottom && left)
	if (chainup (gx, gy, r, gchain, &pcount, b->ch[bcol-1][brow-1]) != 0)
	  return GGLERR_NOMEM;
      if (check_quadrant (gx, gy, gchain, pcount, quad, minpt)) {
	double result = circle_ipol (gx, gy, gchain, pcount, skipg, exponent, eps);
	fprintf (fpgrd, "%12.3f %12.3f %12.3f", gx, gy, result);
	if (dodiff != FALSE) fprintf (fpgrd, " %12.3f %12.3f", gh, gh-result);
	if (entireline[readcount]!='\0') {
	  fprintf (fpgrd, "%s", entireline+readcount);
	} else {
	  fprintf (fpgrd, "\n");
	}
      } else { nancount++; }
      chain_free (gchain, FALSE, pcount);
    } /* if in region */
  } /* while read one line */
  *totpt = totcount;
  return -nancount;
}

/***************************************************************************/
int xyz2grd (FILE *fpxyz, FILE *fppts, const char *grdfn,
	     const double xmin, const double xmax,
	     const double ymin, const double ymax,
	     const double xinc, const double yinc,
	     const double radius, const int exponent, const double eps,
	     const int xn, const int yn,
	     const int surfer, const int quadrant,
	     const int minpt, const int dodiff, const int skipg,
	     const int schlauch, const int verbose)
{
  extern char *program_invocation_short_name;
  struct blockchains *blocks;
  double inpoint[3];
  int col, row;
  long int n=0, total=0;
  struct Grid *rgrd = NULL;
  char entireline[1024];
  clock_t start_time, end_time;
  double elapsed;
  unsigned long ptcount=0;

  start_time = clock();
  /* initialize blockchains and all boundaries */
  blocks = blockchain_init (xmin, xmax, ymin, ymax, xn, yn);
  if (blocks == NULL) {
    fprintf (stderr, "%s: could not allocate memory for blocks\n",
	     program_invocation_short_name);
    return GGLERR_NOMEM;
  }
  
  if (blocks->bxinc < radius || blocks->byinc < radius) {
    fprintf (stderr, "%s: radius smaller than blocksize in x or y\n",
	     program_invocation_short_name);
    blockchain_free (blocks);
    return GGLERR_INVAL;
  }
  if (verbose) { 
    blockchain_info (blocks, FALSE);
    fprintf (stderr, "Radius: %12.3f m, exponent: %d, eps: %g\n",
	     radius, exponent, eps);
    if (quadrant != FALSE) {
      fprintf (stderr, "Using quadrant method.\n");
    } else {
      fprintf (stderr, "Using min. number of points method, min=%d\n", minpt);
    }
    fprintf (stderr, "reading input and filling chains...");
    fflush (stderr);
  }

  /* read input file and fill blockchains */
  while (read_line_static (fpxyz, entireline, 1024) != EOF) {
    if (entireline[0] == '#' || entireline[0] == '$' ||
	entireline[0] == '\n'|| entireline[0] == '\r')
      continue;
    if (sscanf (entireline, " %lf %lf %lf",
		&inpoint[0], &inpoint[1], &inpoint[2]) < 3) {
      blockchain_free (blocks);
      return GGLERR_INOUT;
    }
    if (blockpos(blocks, inpoint[0], inpoint[1], &col, &row) == 0) {
      if (blockchain_insert (blocks, inpoint, col, row)!= 0) {
	fprintf (stderr, "\n%s: could not allocate memory for elements\n",
		 program_invocation_short_name);
	blockchain_free (blocks);
	return GGLERR_NOMEM;
      }
      ptcount++;
    }
  }

  /* initialize output grid */
  if (fppts == NULL) {
    rgrd = grdinit (xmin, xmax, ymin, ymax, xinc, yinc);
    if (rgrd == NULL) {
      fprintf (stderr, "\n%s: could not allocate memory for grid\n",
	       program_invocation_short_name);
      blockchain_free (blocks);
      return GGLERR_NOMEM;
    }
  }
  if (verbose) {
    fprintf (stderr, "done\nfilling grid...");
    fflush (stderr);
  }

  /* loop over output grid and fill it */
  if (fppts == NULL) {
    if ((n=blockchain_fillgrid (blocks,rgrd,radius,exponent,eps,quadrant,minpt))>0) {
      fprintf (stderr, "\n%s: could not allocate memory while filling\n",
	       program_invocation_short_name);
      blockchain_free (blocks);
      return GGLERR_NOMEM;
    }
  } else {
    /* no grid, put input file with target points wanted */
    if ((n=blockchain_procpts (xmin,xmax,ymin,ymax,
			       blocks,fppts,grdfn,radius,exponent,eps,
			       quadrant,minpt,dodiff, skipg, &total)  )>0) {
      switch (n) {
      case GGLERR_NOMEM:
	fprintf (stderr, "\n%s: could not allocate memory while filling\n",
		 program_invocation_short_name);
	blockchain_free (blocks);
	return GGLERR_NOMEM;
      case GGLERR_INOUT:
	fprintf (stderr, "\n%s: I/O error while reading external point file\n",
		 program_invocation_short_name);
	blockchain_free (blocks);
	return GGLERR_INOUT;
      }
      fprintf (stderr, "\n%s: invalid input parameters\n",
	       program_invocation_short_name);
      blockchain_free (blocks);
      return GGLERR_INVAL;
    }
  }

  /* free blockchains */
  blockchain_free (blocks);
  end_time = clock();
  if (verbose) {
    elapsed = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    fprintf (stderr, "done\n%ld points gridded in %.1f seconds.\n",
	     ptcount, elapsed);
    if (fppts == NULL)
      fprintf (stderr, "%ld of %ld points set to NaN.\n",
	       -n, rgrd->nx*rgrd->ny);
    else
      fprintf (stderr, "%ld of %ld points set to NaN.\n", -n, total);
  }

  /* save output grid */
  if (fppts == NULL) {
    g_minmax (rgrd);
    if (surfer!=FALSE) {
      surf_writegrid(grdfn, rgrd, "%10.3f");
    } else {
      if (schlauch != FALSE) {
	g_writeasciigrid (grdfn, rgrd, "%10.3f %10.3f %10,3f\n");
      } else {
	char comment[1024];
	sprintf (comment, "\"%s -r%03.3f -e%05.5f -n%d (%s, %d)\"", 
		 program_invocation_short_name, radius, eps, exponent,
		 (quadrant)?"quadrant":"minpt", (quadrant)?4:minpt);
	gmt_writegrid(grdfn, rgrd, "m", comment, verbose);
      }
    }

  /* free grid */
    g_free (rgrd);
  }

  return EXIT_SUCCESS;
}


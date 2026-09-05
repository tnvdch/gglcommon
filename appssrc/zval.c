/*
 *      FILE: zval.c
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: change time format 
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
 *   AUTHORS: Cocard Marc   <cocard@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            04-06-1999  --  created
 *            30-01-2002  --  major revision and extensions 
 *                            options added etc
 */

#include <errno.h>      /*  the allround headers */
#include <stdlib.h>     /*  the allround headers */
#include <stdio.h>      /*  crucial for file input and output   */
#include <unistd.h>     /*  to parse the command line */
#include <ctype.h>      /*  character handling */
#include <string.h>     /*  string handling */
#include <math.h>       /*  mathematics */

#include "common.h"     /*  our own institute's headers */

#define COPYRIGHT "(c) 2002 Geodesy and Geodynamics Lab, ETH Zurich"

struct Param { 
  int verbose;
  int  rasternr;
  double dummy;
  int prefix_char;
  int prefix_col;
  int switchorder;
  double scalefactor;
  int decimalplaces;
  int interpolationmethod; /* not yet implemented */
};

/* --------------------------------------------------------------------------*/
void display_version(char *pgnname)
{
  fprintf (stderr,"\n%s version %s\n  %s\n", pgnname, SRC_VERSION, COPYRIGHT);
  fprintf (stderr,"  maintainers: M. Cocard <cocard@geod.ethz.ch>\n");
  fprintf (stderr,"\n");
}

/* --------------------------------------------------------------------------
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void display_usage (const int detailed, char *pgnname)
{
  fprintf (stderr, "usage:\n  %s [-hVv] [-i <gridnr>] [-d dummyvalue] [-p <prefixchars>]\n", pgnname);
  fprintf (stderr, "       [-N <pefixcols>] [-f <scalefactor>][ -r <decimalplaces>]\n");
  fprintf (stderr, "\n");
  
  fprintf (stderr, "purpose:\n");
  fprintf (stderr, "  interpolation of rasterfile at given location\n");
  fprintf (stderr, "  reads from stdin -- writes to stdout\n");
  fprintf (stderr, "\n");
  
  if (detailed)
    {
      fprintf (stderr, "options:\n");
      fprintf (stderr,"  -v    verbose\n");
      fprintf (stderr,"  -V    Version string\n");
      fprintf (stderr,"  -h    This help message\n");
      fprintf (stderr,"  -i    nr of rasterfile; see $RASTERINFO (default = 101)\n"); 
      fprintf (stderr,"  -d    dummy value if outside grid (default=-999.9)\n");
      fprintf (stderr,"  -S    switch order of input x/y values\n");
      fprintf (stderr,"  -p    Skip <prefix_char> before parsing values (default = 0)\n");
      fprintf (stderr,"  -N    Skip <prefix_col> before parsing values (default = 0)\n");
      fprintf (stderr,"  -f    scale factor of limits in $RASTERINFO\n");
      fprintf (stderr,"  -r    decimal places of interpolated value\n");
    }  
  else
    fprintf (stderr,"  -h    Prints more detailed help on each option\n\n");
  
  fprintf (stderr, "\n");
}

/* --------------------------------------------------------------------------
 * get information from command lines
 */
int process_command_line(int argc, char **argv, struct Param *param, char *pgnname)
{
  /* argument character to process */
  int arg_char;         
  
  /* Now we process all command line options   */
  while ((arg_char = getopt (argc, argv, "hVvSi:d:p:N:r:f:")) != -1)
    switch (arg_char)
      {
      case 'h': display_usage(TRUE, pgnname);  return -1;              /* -h the help option */
      case 'V': display_version(pgnname);  return -1;                  /* -V the version string */
      case 'v': param->verbose = TRUE; break;
      case 'f': param->scalefactor=strtod(optarg, NULL); break;        /* -g swiss projection */
      case 'd': param->dummy = strtod(optarg, NULL);  break;           /* -d dummy value */   
      case 'p': param->prefix_char = strtol(optarg, NULL, 10); break;  /* -p prefix in chars */
      case 'N': param->prefix_col = strtol(optarg, NULL, 10); break;   /* -N prefix in cols */
      case 'i': param->rasternr = strtol(optarg, NULL, 10); break;     /* -i rasternummer */
      case 'r': param->decimalplaces = strtol(optarg, NULL, 10); break;/* -r decimal places */
      case 'S': param->switchorder = TRUE; break;                      /* -S switching order of xy values */
      case '?': display_usage(FALSE, pgnname); 	return -1;             /* unknown option */
        
      } /* end processing command line option  */
  return 0;
}

/* --------------------------------------------------------------------------*/
int main (int argc, char **argv)
{
  
  struct Grid *g;
  char entireline[MAXLINELENGTH];
  char *pline;
  char format[20];
  double x, y, z, doubledummy;
  long int linecount=0;
  char * pgnname = get_fname_wopath(argv[0]);
  struct Param param;                /* contains all parameters */
  long int datacount = 0, dummycount = 0;

  memset(&param, 0, sizeof(struct Param));
  param.verbose=FALSE;
  param.switchorder=FALSE;
  param.dummy = -999.99;  
  param.rasternr = 101;
  param.scalefactor = -1;
  param.decimalplaces = 2;

  if (process_command_line(argc, argv, &param, pgnname) != 0) return -1;

  
  /* read info concerning rasterfile */
  if ((g = gmt_rasterinfo(param.rasternr)) == NULL)  return -1;
    
  /* GMT caveat in rasterinfofile where limits are in units of 10 km for swiss coordinates */ 
  /* if no scalefactor was given via option we suppose that swiss grids have */
  /*  a rasternr bigger than 100 */
  if (param.scalefactor < 0)        
    param.scalefactor = (param.rasternr < 100)?1:10000.0; 
  g_scalelimits(g, param.scalefactor, param.scalefactor);
  
  /* open rasterfile without reading data */
  if (gmt_openrasterfile(g) != 0) return -1;
 
  if (param.verbose)   
    {
      fprintf(stderr,"\n");
      g_printinfo(stderr, g);
    }

  /* generate format string */
  sprintf(format," %%%d.%df", param.decimalplaces+6, param.decimalplaces);

  while (read_line_static(stdin, entireline, MAXLINELENGTH) != EOF)          /* loop over file */
    {
      linecount++;
      if (iscomment(entireline))   /* copy comments line to output without any changes */
        {
          fprintf (stdout, "%s", entireline);
          continue;
        }

      entireline[strlen(entireline)-1]='\0';
      pline=entireline;
      pline += param.prefix_char;                    /* position s after prefix_char */      
      pline += col_offset(pline, param.prefix_col);  /* position s after prefix_col */
      
      /* get x and y values from string */
      if (sscanf(pline,"%lf %lf", &x, &y) != 2)
	{
	  fprintf(stderr,"error in line %ld\n", linecount);
	  return -1;
	}
      /* exchange x and y if necessary */
      if (param.switchorder)
	{
	  doubledummy = x;
	  x = y;
	  y = doubledummy;
	}

      /* get z-value from raster / dummy if outside of area */
      datacount++;      
      if (g_inter_simple (g, x, y, 4, &z) != 0) 
	{
	  z = param.dummy;
	  dummycount++;
	}

      fprintf(stdout,"%s",entireline);       /* write inputline to output */
      fprintf(stdout,format, z);             /* write z-value to output */
      fprintf(stdout,"\n");
    }

  if (param.verbose)
    {
      fprintf(stderr,"\n");
      fprintf(stderr,"Data points read: %ld\n", datacount);
      fprintf(stderr,"Data points interpolated: %ld (%0.1f%%)\n", 
	      datacount - dummycount,(datacount - dummycount)*100.0/datacount);
      fprintf(stderr,"Data points outside grid: %ld (%0.1f%%)\n", 
	      dummycount, dummycount*100.0/datacount);
      fprintf(stderr,"\n");
    }

  gmt_closerasterfile(g);
  g_free(g);
  return 0;
}

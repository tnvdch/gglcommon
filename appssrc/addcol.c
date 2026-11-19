/*
 *      FILE: addcol
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: merge 2 files using 1 common column
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
 *   AUTHORS: Cocard Marc   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            02-06-1999  --  created (first version)
 * 
 */

#define COPYRIGHT "(c) 1999 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>	/*  to parse the command line */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include "common.h"	/*  our own institute's headers */

#define LINELENGTH 2048
/* maximal length of 1 line */

struct File {
  FILE *fp;
  char *filename;
  long int linecount;
  char entireline[LINELENGTH];
  char *restline;
  double x;
};

/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void display_version()
{
  extern char *program_invocation_short_name; 
  fprintf (stderr,"%s version %s\n  %s\n  %s\n",
	   program_invocation_short_name, SRC_VERSION, LIBVER, COPYRIGHT); 
  fprintf (stderr,"  maintainers: M. Cocard <cocard@geod.ethz.ch>,\n");
  fprintf (stderr,"               E. Favey <favey@geod.ethz.ch>\n\n"); 
  fprintf (stderr," *  DISCLAIMER:   This program is free software; you can redistribute\n");
  fprintf (stderr," *              it and/or modify it under the terms of the GNU General\n");
  fprintf (stderr," *              Public License as published by the Free Software\n");
  fprintf (stderr," *              Foundation; either version 2 of the License, or\n");
  fprintf (stderr," *              (at your option) any later version.\n");
  fprintf (stderr," *                This program is distributed in the hope that it will\n");
  fprintf (stderr," *              be useful, but WITHOUT ANY WARRANTY; without even the\n");
  fprintf (stderr," *              implied warranty of MERCHANTABILITY or FITNESS FOR A\n");
  fprintf (stderr," *              PARTICULAR PURPOSE. See the GNU General Public License\n");
  fprintf (stderr," *              for more details.\n");
  fprintf (stderr," *                You should have received a copy of the GNU General\n");
  fprintf (stderr," *              Public License along with this program; if not, write\n");
  fprintf (stderr," *              to the Free Software Foundation, Inc., 675 Mass Ave,\n");
  fprintf (stderr," *              Cambridge, MA 02139, USA.\n\n");
}										
										       
											
/*
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void display_usage (const int verbose)
{
  extern char *program_invocation_short_name;
  fprintf (stderr,"usage:\n   %s [-hVrdc] [-t tolerance] [-s stringlength] file1 file2\n" , 
	   program_invocation_short_name);
  fprintf (stderr,"\n");
  fprintf (stderr,"  program merges columns from 2 files using the first column of each file\n");
  fprintf (stderr,"  as identifier. (files must always be sorted also when using the -s option).\n");
  fprintf (stderr,"\n\n");
  
  if (verbose != FALSE)
    {
      fprintf (stderr,"  -V    Version string\n");
      fprintf (stderr,"  -h    This help message\n");
      fprintf (stderr,"  -s    identifier taken as string (default identifier double)\n");
      fprintf (stderr,"        you have to provide the length of the string\n");
      fprintf (stderr,"  -t    tolerance [default = 1.0e-6]\n");
      fprintf (stderr,"  -r    print not matching lines to stderr\n");
      fprintf (stderr,"  -c    include comments from first file to output\n");
      fprintf (stderr,"  -d    dummy string for missing data of second file \n");
      fprintf (stderr,"        (in this case all lines from first file are preserved)\n");
      fprintf (stderr,"\n");
    }
}

/*
 * read one line from a file
 */
int readdata(struct File *data, int wrcomment, int stringlength)
{
  do
    {   
      read_line_static(data->fp, data->entireline, LINELENGTH);      
      data->entireline[strlen(data->entireline)-1]='\0';
      data->linecount++;

      if (data->entireline[0]=='#')	/* comment line */
	{
	  if  (wrcomment)
	    fprintf(stdout,"%s\n",data->entireline);
	}
      else	 /* data line */
	{  
	  if (stringlength==0)
	    data->x = strtod(data->entireline,&(data->restline));	  
	  else
	    data->restline = data->entireline+stringlength;
	}
    }
  while (data->entireline[0]=='#' && !feof(data->fp));

  return 0;
}

/*
 * compare two identifiers
 */
int compare(struct File *data1, struct File *data2, double tolerance, int stringlength)
{
  if (stringlength == 0)
    {
      if (fabs(data1->x - data2->x)<tolerance) return 0;
      return (data1->x > data2->x)?1:-1;
    }
  else    
    return strncmp(data1->entireline, data2->entireline, stringlength);
}

/* 
 * get double from option string
*/
int getdouble(char *s, double *x, char *errorstring)
{
  char *ptmp = NULL;
  *x = strtod(s, &ptmp);
  if (ptmp==s)
    {
      fprintf (stderr,"%s:%s\n", program_invocation_short_name, errorstring);
      return -1;
    }
  return 0;
}

/* 
 * get integer from option string
*/
int getint(char *s, int *x, char *errorstring)
{
  char *ptmp = NULL;
  *x = strtol(s, &ptmp, 10);
  if (ptmp==s)
    {
      fprintf (stderr,"%s:%s\n", program_invocation_short_name, errorstring);
      return -1;
    }
  return 0;
}
/*
 * =================================================================
 *  That's where we do the work 
 * ================================================================= 
 */

int main ( int argc, char **argv)
{

  /* argument character to process */
  int arg_char;		

  /* for use in command line parsing */
  double tolerance =1.0e-6;
  int printnonmatching=FALSE;
  int includecomment=FALSE;
  int idstringlength=0;
  char *dummy = NULL;
  struct File data1, data2;
  int order;

  /* for error handling */
  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */
  extern int opterr;                           /* from unistd.h: let's do our own cursing */  

  errinit (argv[0]);
  opterr = 0;

  /*
   * Now we process all command line options
   */
  while ((arg_char = getopt (argc, argv, "hVcd:rt:s:")) != -1)
    switch (arg_char)
      {
      case 'h': display_usage(TRUE);	   /* -h   the help option */
	return EXIT_SUCCESS;

      case 'V': display_version();	   /* -V   the version string */
	return EXIT_SUCCESS;

      case 'd': dummy=strdup(optarg); break;
      case 'r': printnonmatching=TRUE; break;
      case 'c': includecomment=TRUE; break;  
	
      case 't':                             /* -t tolerance for identification */	
	if (getdouble(optarg, &tolerance," Error in tolerance") != 0) return EXIT_FAILURE;	
	break;
 
      case 's':
	if (getint(optarg, &idstringlength, "Error in option -s") != 0) return EXIT_FAILURE;
	break;

	case '?':                              /* unknown option */
	  {
	    if (isprint (optopt))
	      fprintf (stderr, "%s: Unknown option `-%c'.\n",
		       program_invocation_short_name, optopt);
	    else
	      fprintf (stderr,
		       "%s: Unknown option character `\\x%x'.\n",
		       program_invocation_short_name, optopt);
	    display_usage(FALSE);
	    return EXIT_FAILURE;
	  }  
	  
      default: abort (); /* should never be the case */

      } /* end processing command line option  */

  /* if errors in argument list or no argument supported */
  if (argc - optind != 2 || argc < 2) {
    display_usage(FALSE);
    return EXIT_FAILURE;
  }

  memset(&(data1), 0, sizeof(struct File));
  memset(&(data2), 0, sizeof(struct File));

  data1.filename=strdup(argv[optind]);
  data2.filename=strdup(argv[optind+1]);
  
  if ((data1.fp=open_text_file(data1.filename, "r")) == NULL)
  {
	  free(dummy);
	  free( data1.filename);
	  free( data2.filename);
	  return EXIT_FAILURE;
  }
  if ((data2.fp=open_text_file(data2.filename, "r")) == NULL)
  {
 	  free(dummy);
	  free( data1.filename);
	  free( data2.filename);
	  fclose(data1.fp);
	  return EXIT_FAILURE;
  }

  readdata(&data1, includecomment, idstringlength);
  readdata(&data2, FALSE,          idstringlength);
  while (!feof(data1.fp) && !feof(data2.fp))
    {
      if ((order=compare(&data1, &data2, tolerance, idstringlength)) == 0)   /* identifier identical */
	{
	  fprintf(stdout,"%s%s\n",data1.entireline, data2.restline); 
	  readdata(&data1, includecomment,idstringlength);
	  readdata(&data2, FALSE,idstringlength);
	}
      else   /* identifier not equal */
	if (order > 0)  /* get new data from file 2 */
	  {	  
	    if (printnonmatching) fprintf(stderr,"2:%s\n",data2.entireline);
	    readdata(&data2, FALSE, idstringlength);
	  }
	else                            /* get new data from file 1 */
	  {
	    if (dummy != NULL) /* write line filled with dummy string for second file */
	      fprintf(stdout,"%s%s\n",data1.entireline, dummy);
	    else
	      {
		if (printnonmatching) 
		  fprintf(stderr,"1:%s\n", data1.entireline);
	      }
	    readdata(&data1, includecomment,idstringlength);
	  }
    }
  

  free(dummy);
  free( data1.filename);
  free( data2.filename);
  fclose(data1.fp);
  fclose(data2.fp);
  
  return EXIT_SUCCESS;    
}

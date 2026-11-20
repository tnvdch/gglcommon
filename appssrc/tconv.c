/*
 *      FILE: tconv
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
 *   AUTHORS: Cocard Marc   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            20-01-1998  --  created
 * 
 *            02-07-1998  --  line read in static string [1024]
 *                            decomposition of prefix, timestring, postfix 
 *                            simplified (MC)
 *
 *            04-08-1998  --  new features added (MC)
 *                            option -o  time offset
 *                            option -d  number of digits to output
 *
 *            22-14-2000  --  option -N  added (YP)
 *
 *            28-03-2000  --  col_offset used in options -N moved to inout
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>	/*  to parse the command line */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include "common.h"	/*  our own institute's headers */


#define MAXTYP (8)
/* typ = 0     YYYY MM DD HH MM SS.SSS
 * typ = 1     MJD
 * typ = 2     MJD(int) SecOfDay
 * typ = 3     YYYY DayYear SecOfDay
 * typ = 4     WEEK SecofWeek 
 * typ = 5     YYYY DayYear decimalHour
 * typ = 6     JD
 * typ = 7     YYYY.YYY
 * typ = 8     T (Julian centuries since J2000.0)
 */

/* ---------------------------------------------------------------------------------------
 *  Print version information to stderr
 *  You should not need to modify this
 */
void display_version()
{
  extern char *common_version;
  extern char *program_invocation_short_name;
  fprintf (stderr,"%s version %s\n  %s\n  %s\n",
	   program_invocation_short_name, SRC_VERSION, common_version, COPYRIGHT);
  fprintf (stderr,"  maintainers: M. Cocard <cocard@geod.ethz.ch>,");
  fprintf (stderr," E. Favey <favey@geod.ethz.ch>\n\n");
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

/* ---------------------------------------------------------------------------------------
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void display_usage (const int verbose)
{
  extern char *program_invocation_short_name;
  fprintf (stderr,"usage:\n   %s [-hVvk] [-p <prefix_char>] [-N <prefix_col>]\n" , 
	   program_invocation_short_name);
  fprintf (stderr,"         [-I <input formtype>] [-O <output formtype>]\n");
  fprintf (stderr,"         [-o <timeoffset(sec)>] [-d #digits]\n\n");
  fprintf (stderr,"   conversion of different time representations\n");
  fprintf (stderr,"   reads from standard in and writes to standard out\n\n");
  
  if (verbose != FALSE)
    {
      fprintf (stderr,"  -V    Version string\n");
      fprintf (stderr,"  -h    This help message\n");
      fprintf (stderr,"  -v    verbose mode, log conversion progress to stderr\n");
      fprintf (stderr,"  -k    keep input time and insert output type\n");
      fprintf (stderr,"  -p    Skip <prefix_char> before parsing values (default = 0)\n");
      fprintf (stderr,"  -N    Skip <prefix_col> before parsing values (default = 0)\n");
      fprintf (stderr,"  -I    formtype on input\n");
      fprintf (stderr,"  -O    formtype on output\n");
      fprintf (stderr,"  -o    offset in seconds to be added to the output time (default = 0)\n");
      fprintf (stderr,"  -d    number of digits to be printed (default = 6 for seconds, 9 for days\n");
      fprintf (stderr,"        and 15 for centuries)\n");
      fprintf (stderr,"\n");
    }
  else
    fprintf (stderr,"  -h    Prints more detailed help on each option\n\n");      

  fprintf (stderr, "supported formtypes:\n");
  fprintf (stderr, "   0    year month day hour min sec [default]\n");
  fprintf (stderr, "   1    mjd(double)\n");
  fprintf (stderr, "   2    mjd(int) SecOfDay\n");
  fprintf (stderr, "   3    year DayOfYear SecOfDay\n");
  fprintf (stderr, "   4    GPSWeek SecOfWeek\n");
  fprintf (stderr, "   5    year DayOfYear HourOfDay\n");
  fprintf (stderr, "   6    julian date\n");
  fprintf (stderr, "   7    year(double)\n");
  fprintf (stderr, "   8    T(double) (Julian centuries since J2000.0)\n");
  fprintf (stderr,"\n");
}

/* ---------------------------------------------------------------------------------------  
 *  time conversion of string to modified julian date mjd 
 */
int convertstring2time( const  char * tstr, const int inptyp, double * mjd, int *count)
{
  unsigned int mjd0, week, dayofyear;
  double secofday, secofweek, year, fracyear, T;
  struct ttime t;

  switch (inptyp)
    {
    case 0: /* YYYY MM DD HH MM SS.SSSS */
      {
	if (sscanf(tstr,"%d %d %d %d %d %lf %n", 
		   &(t.year), &(t.month), &(t.day), 
		   &(t.hour), &(t.min), &(t.sec), count) != 6) return -1;
	datmjd(&t);
	*mjd =t.mjd;
      }
      break;

    case 1: /* MJD (double) */
      {
	if (sscanf(tstr,"%lf %n", mjd, count) != 1) return -1;
      }
      break;

    case 2: /* MJD (int) SecOfDay */
      {
 	if (sscanf(tstr,"%d %lf %n", &mjd0, &secofday, count) != 2) return -1;
	*mjd = (double) mjd0 + (secofday/86400.0);
      }
     break;

    case 3: /*  YYYY DayofYear SecOfDay */
      {
 	if (sscanf(tstr,"%d %d %lf %n", &(t.year), &dayofyear, &secofday, count) != 3) return -1;
	*mjd = yearday_to_mjd (t.year, dayofyear);
        *mjd += secofday/86400.0;
      }
      break;

    case 4: /*  Week SecofWeek */
      {
	if (sscanf(tstr,"%d %lf %n",  &week, &secofweek, count) != 2) return -1;
	*mjd = gps_to_mjd (secofweek, week);
      }
      break;

    case 5: /*  YYYY DayofYear decimalHour */
      {
	      double hourofday;
 	if (sscanf(tstr,"%d %d %lf %n", &(t.year), &dayofyear, &hourofday, count) != 3) return -1;
	*mjd = yearday_to_mjd (t.year, dayofyear);
        *mjd += hourofday/24.0;
      }
      break;

    case 6: /* JD (double) */
      {
	if (sscanf(tstr,"%lf %n", mjd, count) != 1) return -1;
	*mjd -= JULDAT0;
      }
      break;

    case 7: /* year(double) */
      {
	if (sscanf(tstr,"%lf %n", &year, count) != 1) return -1;
	fracyear=modf(year, &year);
	t.year = (int) year;
	*mjd = yearday_to_mjd (t.year, 1);
	*mjd += fracyear*((leapyear(t.year))?366.0:365.0);
      }
    case 8: /* T(double) */
    {
        if (sscanf(tstr, "%lf %n", &T, count) != 1) return -1;
        *mjd = T * JULCENTURY + MJD_J2000;
    }


   } /* end switch */

  /* reincludes leading blanks / tstr[count] points to start of postfix */ 
  while (*count>0 && (tstr[*count-1]==' ' || tstr[*count-1]=='\n')) (*count)--;

  return 0;
}

/* ---------------------------------------------------------------------------------------   
 * write time to file
 */
int writetime( const int outtyp, double mjd, int digits)
{
  int year, day, week;
  double seconds, yeardouble, dummy;
  int digit_year, digit_hour, digit_sec, digit_cent;
  char format[MAXLINELENGTH];


  digit_year = (digits<0)?9:digits;
  digit_hour = (digits<0)?5:digits;
  digit_sec  = (digits<0)?6:digits;
  digit_cent = (digits<0)?15:digits;

  switch (outtyp)
    {
    case 0: /* YYYY MM DD HH MM SS.SSSS */
      {
	char *ptmp = strtime2(mjd, digit_sec);
	fprintf(stdout,"%s", ptmp);
	free (ptmp);
      }
      break;
      
    case 1: /* MJD (double) */
      {
	sprintf(format,"%%%d.%df",digit_year+6, digit_year);
	fprintf(stdout,format, mjd);
      }
      break;

    case 2: /* MJD (int) SecOfDay */
      {
	seconds = (mjd-floor(mjd))*86400.0;
	sprintf(format,"%%5.0f %%%d.%df",digit_sec+6, digit_sec);
	fprintf(stdout, format, floor(mjd), seconds);
      }
      break;
      
    case 3: /*  YYYY DayofYear SecOfDay */
      {
	mjd_to_yearday (mjd, &year, &day);
	seconds = (mjd-floor(mjd))*86400.0;
	sprintf(format,"%%4d %%3d %%%d.%df",digit_sec+5, digit_sec);
	fprintf(stdout,format, year, day, seconds);            
      }
      break;

    case 4: /*  Week SecofWeek */
      {
	mjd_to_gps(mjd, &seconds, &week);
	sprintf(format,"%%4d %%%d.%df",digit_sec+7, digit_sec);
	fprintf(stdout,format, week, seconds);
      }
      break;
      
    case 5: /*  YYYY DayofYear decimalHour */
    {
      double hour;
      mjd_to_yearday (mjd, &year, &day);
      hour = (mjd-floor(mjd))*24.0;
      sprintf(format,"%%4d %%3d %%%d.%df",digit_hour+3, digit_hour);
      fprintf(stdout,format, year, day, hour);
    }
    break;
    
    case 6: /* JD (double) */
      {
	sprintf(format,"%%%d.%df",digit_year+6, digit_year);
	fprintf(stdout,format, mjd+JULDAT0);
      }
      break;
      
    case 7: /* year (double) */
      {
	mjd_to_yearday(mjd, &year, &day);
	yeardouble=year + (day-1+modf(mjd, &dummy))/(leapyear(year)?366.0:365.0);
	sprintf(format,"%%%d.%df",digit_year+5, digit_year);
	fprintf(stdout,format,yeardouble);
      }
      break;

    case 8: /* T(double) */
    {
        sprintf(format, "%%.%df", digit_cent);
        fprintf(stdout, format, (mjd-MJD_J2000)/JULCENTURY);
    }

      
    } /* end switch */
  
  return 0;
}

/* --------------------------------------------------------------------------------------- 
 *  That's where we do the work 
 */

int main ( int argc, char **argv)
{

  /* argument character to process */
  int arg_char;		

  /* default initialization of param used in command line parsing */
 // int verbose = FALSE;
  int keep = FALSE;
  int prefix_char = 0;
  int prefix_col = 0;
  int inptyp = 0;
  int outtyp = 0;
  double timeoffset = 0.0;
  int digits = -1;

 /* for use in main */
  char *ptmp;
  char entireline[MAXLINELENGTH];
  char *s;  
  double mjd;
  unsigned long linecount = 0;
  int charcount;
  int offset; /* Char Offset resulting from the overreadinf of t_column */


  /* for error handling */
  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */
  extern int opterr;                           /* from unistd.h: let's do our own cursing */  

  errinit (argv[0]);
  opterr = 0;
  
  /*
   * Now we process all command line options
   */
  while ((arg_char = getopt (argc, argv, "hVvkp:I:O:o:d:N:")) != -1)
    switch (arg_char)
      {
      case 'h': display_usage(TRUE);	   /* -h   the help option */
	      return EXIT_SUCCESS;

      case 'V': display_version();	   /* -V   the version string */
	return EXIT_SUCCESS;

//      case 'v': verbose = TRUE;       	break;
	
      case 'k': keep = TRUE;       	break;
	
      case 'o': timeoffset = strtod(optarg, &ptmp); 	break;  /* -o   timeoffset */

      case 'd': digits = strtol(optarg, &ptmp, 10); break;    /*  -d digits */

      case 'p':                            /* -p  length of prefix */
	{
	  char *ptmp = NULL;
	  prefix_char = strtol(optarg, &ptmp, 10);
	  if (ptmp == optarg || prefix_char < 0)
	    {
	      fprintf (stderr,"%s: Option `-p' requires positive integer argument.\n", 
		       program_invocation_short_name);
	      display_usage(FALSE);
	      return EXIT_FAILURE;
	    }
	}
	break;

      case 'N':                            /* -N # of columns to skip */
	{
	  char *ptmp = NULL;
	  prefix_col = strtol(optarg, &ptmp, 10);
	  if (ptmp == optarg || prefix_col < 0)
	    {
	      fprintf (stderr,"%s: Option `-N' requires positive integer argument.\n", 
		       program_invocation_short_name);
	      display_usage(FALSE);
	      return EXIT_FAILURE;
	    }
	}
	break;
 	
      case 'I':  /* -I # of input format of time string */
	{
	  char *ptmp = NULL;
	  inptyp = strtol(optarg, &ptmp, 10);
	  if (ptmp == optarg || inptyp < 0 || inptyp > MAXTYP)
	    {
	      fprintf (stderr,"%s: Option `-I' requires positive integer argument in the range [0 - %d]\n", 
		       program_invocation_short_name, MAXTYP);
	      display_usage(FALSE);
	      return EXIT_FAILURE;
	    }
	} 
	break;
	
      case 'O': /* -O  # of output format of time string */
	{
	  outtyp = strtol(optarg, &ptmp, 10);
	  if (ptmp == optarg || outtyp < 0 || outtyp > MAXTYP)
	    {
	      fprintf (stderr,"%s: Option `-O' requires positive integer argument in the range [0 - %d]\n", 
		       program_invocation_short_name, MAXTYP);
	      display_usage(FALSE);
	      return EXIT_FAILURE;
	    }
	}
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
  if (argc - optind != 0 || argc == 1) {
    display_usage(FALSE);
    return EXIT_FAILURE;
  }

  while (read_line_static(stdin, entireline, MAXLINELENGTH) != EOF) {         /* loop over file */
    linecount++;
    if (iscomment(entireline))   /* copy comments line to output without any changes */
	{
	  fprintf (stdout, "%s", entireline);
	  continue;
	}
	  
    if (prefix_char >= strlen(entireline)) { /* error in inputline */
	fprintf (stderr, "%s: invalid input in line %ld (prefix longer than line)\n",
		 program_invocation_short_name, linecount);
	return EXIT_FAILURE;
    }

    s=entireline + prefix_char; /* position s after prefix_char */

    if ((offset=col_offset(s, prefix_col)) < 0) return  EXIT_FAILURE; /* position s after prefix_col */
    s += offset;

    if (convertstring2time(s, inptyp, &mjd, &charcount) != 0) { /* error in conversion */
	fprintf (stderr, "%s:error in converting line nr %ld\n",
		 program_invocation_short_name, linecount);
	return EXIT_FAILURE;
    }
      
    entireline[prefix_char + offset]='\0';  /* clip entire line to prefix */
    fprintf(stdout,"%s ",entireline);  /* write prefix */
    mjd+=(timeoffset/86400);           /* add time offset */
    writetime(outtyp, mjd, digits);    /* write time info */
    if (keep) {
	    fprintf(stdout," ");
	    writetime(inptyp, mjd, digits);  /* keep old time info */
    }
    s+=charcount;                    /* position to start of postfix (always contains at least \n */
    fprintf(stdout,"%s",s);          /* write postfix */

  }
  return EXIT_SUCCESS;    
}

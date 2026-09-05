/*
 *      FILE: rayp
 *
 *   PROJECT: Laser Altimetry on Glacier
 *
 *   PURPOSE: calculate scanning width
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            14-07-1997  --  created
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1997 Geodesy and Geodynamics Lab, ETH Zurich"

#define INITSTATE (0)
#define CALCSTATE (1)

/*
 *  We need some standard headers:
 *
 *  the allround headers
 */
#include <errno.h>
#include <stdlib.h>
/*
 *  this one is crucial for file input and output  
 */
#include <stdio.h>
/*
 *  this one to parse the command line
 */
#include <unistd.h>
/*
 *  these to do character and string handling
 */
#include <ctype.h>
#include <string.h>
/*
 *  provide some error handling that is GNU extension
 *  (should be, I do have an older system, so I had to hack it)
 */
#include <math.h>
/*
 *  our own institute's headers
 */
#include <common.h>

/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void
display_version()
{
   extern char *common_version;
   extern char *program_invocation_short_name;
  fprintf (stderr,
	   "%s version %s \n  %s\n  %s\n",
	   program_invocation_short_name, SRC_VERSION, common_version,
	   COPYRIGHT);
       fprintf (stderr,
       "  maintainer: E. Favey <favey@geod.ethz.ch>\n\n");
}

/*
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void
display_usage (
	       const int verbose
	       )
{
   extern char *program_invocation_short_name;
   fprintf (stderr,
   "usage:\n   %s [-Vhvmt] [-a <angle>] [-f <bln-flag>] [-u <altitude>]\n     input_file output_file\n\n" ,
	    program_invocation_short_name);
   if (verbose == 0)
     {
       fprintf (stderr,
		"  -V    Version string\n");
       fprintf (stderr,
		"  -h    This help message\n");
       fprintf (stderr,
		"  -v    verbose mode, tell about conversion progress\n");
       fprintf (stderr,
		"  -f    BLN flag to use on output polygons, default = 0\n");
       fprintf (stderr,
		"  -t    multiple files mode, one trajectory per file\n");
       fprintf (stderr,
		"  -m    multiple files mode, one polygon per file\n");
       fprintf (stderr,
		"  -a    scanning angle [degrees]\n");
       fprintf (stderr,
		"  -u    altitude, when not provided in input file\n");
       fprintf (stderr,
       "  <file>   Files read and written, should be BLN format\n");
     }
   else
     {
        fprintf (stderr,
		"  -h    Prints more detailed help on each option\n\n");      
     }
}

int scanwidth (
  FILE *fpin,
  const char *outname,
  const double scanangle, /* in rad */
  const double altitude,
  const int verbose,
  const int writemode,
  const int outflag)
{
  int i;
  char filename[256];
  FILE *fpout;
  /* two vectors pointing to trajectory */
  double p1[3], p2[3];
  /* tangent vector along flight path */
  double fp[3];
  /* vector along local vertical */
  double lv[3], olv[3];
  /* vector to ground point below local vertical */
  double gp[3], ogp[3];
  /* vectors to the target points of laser beam */
  double t11[3], t12[3], t21[3], t22[3];
  /* radius vector from ground point to target points */
  double er[3]; /* unit vector in direction of r */
  double mr; /* magnitude of r */
  double r[3];
  /* altitude of craft */
  double aoc;
  /* altitude of ground (could be read from separate grid) */
  double aog;
  /* number of coordinate pairs that follow */
  int nocp;
  /* flag in BLN file */
  int blnflag;
  /* string to hold entire input line */
  char *theline;
  /* read state */
  int state = INITSTATE;
  int cols;
  int rows;

  if (writemode == 0)
    {
      if (verbose) fprintf (stderr, "%s: creating %s ...\n",
			    program_invocation_short_name, outname);
      if ( !(fpout = fopen (outname, "w")) )
	{
	  fprintf (stderr,
		   "%s: Could not open file `%s' for writing: %s\n",
		   program_invocation_short_name,
		   outname,
		   strerror (errno) );
	  return EOF;
	}
    }
  i = 0;
  while ((read_line (fpin, &theline)) != EOF)
    {
      switch (state)
	{
	case INITSTATE:
	  /* we accept empty lines preceeding header */
	  if (strlen (theline) < 4)
	    {
	      free (theline);
	      continue;
	    }
	  /* read header info: length and flag */
	  sscanf (theline, "%d %d", &nocp, &blnflag);
	  free (theline);
	  rows = 0;
	  /* read first point and treat is as last read */
	  if (read_line (fpin, &theline) == EOF)
	    {
	      fprintf (stderr, "%s: corrupt input file, no BLN spec\n",
		       program_invocation_short_name);
	      fprintf (stderr, "%s: aborting program execution ...\n",
		       program_invocation_short_name);
	      return EOF;
	    }
	  cols = sscanf (theline, "%lf %lf %lf %lf",
			 &(p2[0]), &(p2[1]), &aog, &aoc);
	  rows++;
	  if (cols < 3)
	    {
	      fprintf (stderr, 
		       "%s: input file does not contain third column\n",
		       program_invocation_short_name);
	      fprintf (stderr, 
		       "%s: consult manpage for problem description\n",
		       program_invocation_short_name);
	      fprintf (stderr, "%s: aborting program execution ...\n",
		       program_invocation_short_name);
	      return EOF;
	    }
	  /* read grid info into aog and aoc */
	  /* to be extended */
	  if (altitude)
	    {
	      if (verbose)
		fprintf (stderr,
			 "%s: using altitude provided on commandline\n",
			 program_invocation_short_name);
	      aoc = p2[2] = altitude;
	    }
	  else
	    {
	      if (cols == 4)
		{
		  if (verbose)
		    fprintf (stderr,
		    "%s: using craft altitude provided in input file\n",
			     program_invocation_short_name);
		  p2[2] = aoc;
		}
	      else
		{
		  fprintf (stderr, 
			   "%s: input file does not contain fourth column\n",
			   program_invocation_short_name);
		  fprintf (stderr, 
			   "%s: consult manpage for problem description\n",
			   program_invocation_short_name);
		  fprintf (stderr, "%s: aborting program execution ...\n",
			   program_invocation_short_name);
		  return EOF;
		}
	    }
	  free (theline);
	  /* calculate local vertical newpoint */
	  lv[0] = lv[1] = 0.0;
	  lv[2] = aog - aoc;
	  /* calculate ground point */
	  v3_add (p2, lv, gp);
	  state = CALCSTATE;
	  if (writemode == 2)
	    {
	      sprintf (filename, "%s%03d.bln\0", outname, i++);
	      if (verbose) fprintf (stderr, "%s: creating %s ...\n",
				    program_invocation_short_name, filename);
	      if ( !(fpout = fopen (filename, "w")) )
		{
		  fprintf (stderr,
			   "%s: Could not open file `%s' for writing: %s\n",
			   program_invocation_short_name,
			   filename,
			   strerror (errno) );
		  return EOF;
		}
	    }
	  break;
	case CALCSTATE:
	  /* move last read to initial trajectory point */
	  v3_copy (p2, p1);
	  v3_copy (lv, olv);
	  v3_copy (gp, ogp);
	  /* extract needed values from theline */
	  cols = sscanf (theline, "%lf %lf %lf %lf",
			 &(p2[0]), &(p2[1]), &aog, &aoc);
	  if (cols < 3)
	    {
	      fprintf (stderr, 
		       "%s: input file does not contain third column\n",
		       program_invocation_short_name);
	      fprintf (stderr, 
		       "%s: consult manpage for problem description\n",
		       program_invocation_short_name);
	      fprintf (stderr, "%s: aborting program execution ...\n",
		       program_invocation_short_name);
	      return EOF;
	    }
	  /* read grid info into aog and aoc */
	  /* to be extended */
	  if (altitude)
	    {
	      if (verbose)
		fprintf (stderr,
			 "%s: using altitude provided on commandline\n",
			 program_invocation_short_name);
	      aoc = p2[2] = altitude;
	    }
	  else
	    {
	      if (cols == 4)
		{
		  if (verbose)
		    fprintf (stderr,
		    "%s: using craft altitude provided in input file\n",
			     program_invocation_short_name);
		  p2[2] = aoc;
		}
	      else
		{
		  fprintf (stderr, 
			   "%s: input file does not contain fourth column\n",
			   program_invocation_short_name);
		  fprintf (stderr, 
			   "%s: consult manpage for problem description\n",
			   program_invocation_short_name);
		  fprintf (stderr, "%s: aborting program execution ...\n",
			   program_invocation_short_name);
		  return EOF;
		}
	    }
	  free (theline);
	  /* calculate direction of flight */
	  v3_sub (p2, p1, fp);
	  /* calculate local vertical newpoint */
	  lv[0] = lv[1] = 0.0;
	  lv[2] = aog - aoc;
	  /* calculate ground point */
	  v3_add (p2, lv, gp);
	  /* calculate target points */
	  v3_cross (fp, olv, r);
	  v3_norm (r, er);
	  mr = v3_length (olv) * tan (scanangle);
	  v3_multx (er, mr, r);
	  v3_add (ogp, r, t11);

	  v3_multx (er, -mr, r);
	  v3_add (ogp, r, t12);

	  v3_cross (fp, lv, r);
	  v3_norm (r, er);
	  mr = v3_length (lv) * tan (scanangle);
	  v3_multx (er, -mr, r);
	  v3_add (gp, r, t21);

	  v3_multx (er, mr, r);
	  v3_add (gp, r, t22);
	  /* write closed area with target points */
	  /* fprintf (fpout, "5 1\n");
	   * fprintf (fpout, "%12.3f %12.3f %12.3f\n",
	   *    t11[0], t11[1], t11[2]);
	   * fprintf (fpout, "%12.3f %12.3f %12.3f\n",
	   *    t12[0], t12[1], t12[2]);
	   * fprintf (fpout, "%12.3f %12.3f %12.3f\n",
	   *    t21[0], t21[1], t21[2]);
	   * fprintf (fpout, "%12.3f %12.3f %12.3f\n",
	   *    t22[0], t22[1], t22[2]);
	   * fprintf (fpout, "%12.3f %12.3f %12.3f\n",
	   *    t11[0], t11[1], t11[2]);
	   * fprintf (fpout, "\n");
	   */
	  /* write closed area with target points */
	  if (writemode == 1)
	    {
	      sprintf (filename, "%s%03d.bln\0", outname, i++);
	      if (verbose) fprintf (stderr, "%s: creating %s ...\n",
				    program_invocation_short_name, filename);
	      if ( !(fpout = fopen (filename, "w")) )
		{
		  fprintf (stderr,
			   "%s: Could not open file `%s' for writing: %s\n",
			   program_invocation_short_name,
			   filename,
			   strerror (errno) );
		  return EOF;
		}
	    }
	  fprintf (fpout, "5,%d\n", outflag);
	  fprintf (fpout, "%-.3f,%-.3f\n",
		   t11[0], t11[1], t11[2]);
	  fprintf (fpout, "%-.3f,%-.3f\n",
		   t12[0], t12[1], t12[2]);
	  fprintf (fpout, "%-.3f,%-.3f\n",
		   t21[0], t21[1], t21[2]);
	  fprintf (fpout, "%-.3f,%-.3f\n",
		   t22[0], t22[1], t22[2]);
	  fprintf (fpout, "%-.3f,%-.3f\n",
		   t11[0], t11[1], t11[2]);
	  if (writemode == 1) fclose (fpout);
	  if (++rows >= nocp)
	    {
	      state = INITSTATE;
	      if (writemode == 2) fclose (fpout);
	      if (verbose)
		fprintf (stderr, "%s: end of subset\n",
			 program_invocation_short_name);
	    }
	  /* stay in CALCSTATE */
	default:
	  /* curse badly, made programming mistake */
	}
      /* if length is reached, read new line or quit */
    }
  if (writemode == 0) fclose (fpout);
  return 0;
}

/*
 *  That's where we do the work 
 */
int
main (
      int argc,        /* Number of arguments         */
      char **argv      /* array of argument strings   */
      )
{
  int
    arg_index,         /* index variable for loop       */
    arg_char;          /* argument character to process */
  int use_pipe = FALSE;
  int verbose = FALSE;
  int writemode = 0;
  double scanangle = 13.4; /* INS SCALARS II */
  double altitude = 0.0;
  int outflag = 0;
  char *ptmp;
  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */

  extern int opterr;      /* from unistd.h: let's do our own cursing */

  errinit (argv[0]);      /* from errproc.h */
  opterr = 0;

  /*
   * Now we process all command line options
   */
  while ((arg_char = getopt (argc, argv, "Vhvmta:u:f:")) != -1)
    switch (arg_char)
      {
      case 'h':             /* -h   the help option */
	display_usage(0);
      case 'V':             /* -v   the version string */
	display_version();
	return EXIT_SUCCESS;
      case 'f':             /* -f <bln-flag> */
        ptmp = NULL;
        outflag = strtol(optarg, &ptmp, 10);
        if (ptmp == optarg || outflag < 0)
          {
            fprintf (stderr, 
                     "%s: Option `-f' requires positive integer argument.\n",
                     program_invocation_short_name);
            display_usage(1);
            return EXIT_FAILURE;
          }
        break;
      case 'a':
	ptmp = NULL;
	scanangle = strtod(optarg, &ptmp);
	if (ptmp == optarg || scanangle < 0)
	  {
	    fprintf (stderr, 
		     "%s: Option `-a' requires positive double argument.\n",
		     program_invocation_short_name);
	    display_usage(1);
	    return EXIT_FAILURE;
	  }
	break;
      case 'u':
	ptmp = NULL;
	altitude = strtod(optarg, &ptmp);
	if (ptmp == optarg || altitude < 0)
	  {
	    fprintf (stderr, 
		     "%s: Option `-u' requires positive double argument.\n",
		     program_invocation_short_name);
	    display_usage(1);
	    return EXIT_FAILURE;
	  }
	break;
      case 'v':
	verbose = TRUE;
	break;
      case 'm':
	if (writemode != 0)
	  {
	    fprintf (stderr,
	    "%s: Option `-m' may not be used in conjuntion with `-t'.\n",
		   program_invocation_short_name);
	    return EXIT_FAILURE;
	  }
	writemode = 1;
	break;
      case 't':
	if (writemode != 0)
	  {
	    fprintf (stderr,
	    "%s: Option `-m' may not be used in conjuntion with `-t'.\n",
		   program_invocation_short_name);
	    return EXIT_FAILURE;
	  }
	writemode = 2;
	break;
      case '?':             /* oops, an unknown option */
	if (isprint (optopt))
	  fprintf (stderr, "%s: Unknown option `-%c'.\n",
		   program_invocation_short_name, optopt);
	else
	  fprintf (stderr, "%s: Unknown option character `\\x%x'.\n",
		   program_invocation_short_name, optopt);
	display_usage(1);
	return EXIT_FAILURE;
      default:              /* actually, this should never happen */
	abort ();
      }

  /*
   *  Check if we do have all the necessary items
   */
  if (  ((argc - optind) < 2) && use_pipe == FALSE )
    {
      display_usage(1);
      return EXIT_SUCCESS;     /* If there were no arguments ... */
    }


  /*
   * Now we parse the non-option arguments
   */
    {
      /* Those variables are only needed locally */
      FILE *fpin, *fpout;


      if ( !(fpin = fopen (argv[optind], "r")) )
	{
	  fprintf (stderr,
		   "%s: Could not open file `%s' for reading: %s\n",
		   program_invocation_short_name,
		   argv[optind],
		   strerror (errno) );
	  return EXIT_FAILURE;
	}
/*       if ( !(fpout = fopen (argv[optind+1], "w")) ) */
/* 	{ */
/* 	  fprintf (stderr, */
/* 		   "%s: Could not open file `%s' for writing: %s\n", */
/* 		   program_invocation_short_name, */
/* 		   argv[optind+1], */
/* 		   strerror (errno) ); */
/* 	  return EXIT_FAILURE; */
/* 	} */

      /*
       *  this is where we do all the action!
       *  read all kinds of stuff from file
       */
      scanangle *= M_PI/180.0; /* convert degrees to rad */
      scanwidth (fpin, argv[optind+1], scanangle, altitude,
		 verbose, writemode, outflag);

      /* clean up */
      fclose (fpin);
/*       fclose (fpout); */
    }
  /*
   * We are all done :-)
   * Exit properly!
   */
  return EXIT_SUCCESS;
}


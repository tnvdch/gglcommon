/*
 *      FILE: 
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: Calculate scan path on ground
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            17-05-1997  --  created
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1997 Geodesy and Geodynamics Lab, ETH Zurich"

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
void display_version()
{
	extern char *program_invocation_short_name;
	extern char *common_version;
	fprintf (stderr,
		 "%s version %s\n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION,
		 common_version, COPYRIGHT);
	fprintf (stderr,
		 "  maintainers: E. Favey <favey@geod.ethz.ch>\n\n");
}

/*
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void display_usage (
	const int verbose
	)
{
	extern char *program_invocation_short_name;
	fprintf (stderr,
   "usage:\n   %s [-Vh] [-i <maxiter>] [-a <scan angle>] <file.nav>\n\n" ,
	    program_invocation_short_name);
	if (verbose == TRUE)
	{
       fprintf (stderr,
		"  -V              Display version information\n");
       fprintf (stderr,
		"  -h              This help message\n");
       fprintf (stderr,
		"  -a              Scanning angle adding to roll\n");
       fprintf (stderr,
		"  -i              maximum # of iterations per point\n");
       fprintf (stderr,
		"  <file.nav>      nav file to process, swiss coordinates\n");
	}
	else
	{
        fprintf (stderr,
		"  -h    Prints more detailed help on each option\n\n");      
	}
}







/*
 *  That's where we do the work 
 */
int main (
	int argc,        /* Number of arguments         */
	char **argv      /* array of argument strings   */
	)
{
  int
	  /*  arg_index,          index variable for loop       */
    arg_char;          /* argument character to process */

  int maxiter = 20;
  double sangle = 0.0;

  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */

  extern int opterr;      /* from unistd.h: let's do our own cursing */

  errinit (argv[0]);      /* from errproc.h */
  opterr = 0;

  /*
   * Now we process all command line options
   */
  while ((arg_char = getopt (argc, argv, "Vhi:a:")) != -1)
    switch (arg_char)
      {
      case 'h':             /* -h   the help option */
	display_usage(TRUE);
      case 'V':             /* -v   the version string */
	display_version();
	return EXIT_SUCCESS;

      case 'i':
        {
          char *ptmp = NULL;
          maxiter = strtol(optarg, &ptmp, 10);
          if (ptmp == optarg || maxiter < 0)
            {
              fprintf (stderr, 
                       "%s: Option `-i' requires positive integer argument.\n",
                       program_invocation_short_name);
              display_usage(FALSE);
              return EXIT_FAILURE;
            }
        }
        break;
      case 'a':
        {
          char *ptmp = NULL;
          sangle = strtod(optarg, &ptmp);
          if (ptmp == optarg)
            {
              fprintf (stderr, 
                       "%s: Option `-a' requires double argument.\n",
                       program_invocation_short_name);
              display_usage(FALSE);
              return EXIT_FAILURE;
            }
        }
        break;

      case '?':             /* oops, an unknown option */
	if (isprint (optopt))
	  fprintf (stderr, "%s: Unknown option `-%c'.\n",
		   program_invocation_short_name, optopt);
	else
	  fprintf (stderr, "%s: Unknown option character `\\x%x'.\n",
		   program_invocation_short_name, optopt);
	display_usage(FALSE);
	return EXIT_FAILURE;
      default:              /* actually, this should never happen */
	abort ();
      }

  /*
   *  Check if we do have all the necessary items
   */
  if ( (argc - optind) != 1 )
    {
      display_usage(FALSE);
      return EXIT_SUCCESS;     /* If there were no arguments ... */
    }


  /*
   * Now we parse the non-option arguments
   */
    {
	    FILE *fpnav;
	    struct Grid *dtm;
	    double pos[3], att[3], ground[3], tmp1[3];
	    double distance;
	    double secofday;
	    int count, i;
	    struct Angle tmpa;
		int satcount;

	    tmpa = afromdeg (sangle);
	    sangle = tmpa.rad;

	    dtm = gmt_rasterinfo(100); /* DHM 25 */
	    if (dtm==NULL) return EXIT_FAILURE;
	    if (gmt_openrasterfile (dtm) != 0)
		    return EXIT_FAILURE;
	    g_printinfo (stderr, dtm);

	    fpnav = open_text_file (argv[optind+1], "r");
	    if (fpnav == NULL) return EXIT_FAILURE;
	    
	    fprintf (stdout, "$ secofday  LK Y      LK X     height    distance iter\n");

	    /* Data must be ch1903+, bessel elliptic */
	    while (feof (fpnav) == 0) {
		    count = fscanf (fpnav, " %lf %lf %lf %lf %lf %lf %lf %d",
				    &secofday, &(pos[1]), &(pos[0]), &(pos[2]),
				    &(att[2]), &(att[1]), &(att[0]), &satcount);

		    /* convert all angles to rad */
		    for (i=0;i<3;i++) {
			    tmpa = afromdeg(att[i]);
			    att[i] = tmpa.rad;
		    }
	/*	    for (i=0;i<2;i++) {
			    tmpa = afromdeg(pos[i]);
			    pos[i] = tmpa.rad;
		    }
	*/
	/*	    ch_elltomil(pos[0], pos[1], &(pos[1]), &(pos[0]));
	*/	    
		    att[0] = sangle + att[0];
		    count = g_intersect_vector (dtm, maxiter,
						pos, att, ground);
		    if (count < 0) {
#ifdef WITH_DEBUG
			    fprintf (stderr, "%s: datapoint not above DTM\n",
				     program_invocation_short_name);
#endif
			    continue;
		    }

		    v3_sub (pos, ground, tmp1);
		    distance = v3_length (tmp1);
		    fprintf (stdout, "%12.5f %10.3f %10.3f %10.3f %8.3f %d\n",
			     secofday, ground[0], ground[1], ground[2],
			     distance, count);
	    }
	    fclose (fpgrd);
	    fclose (fpnav);
    }
  /*
   * We are all done :-)
   * Exit properly!
   */
  return EXIT_SUCCESS;
}






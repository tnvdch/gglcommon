/*
 *      FILE: $Source$
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: main program
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Favey Etienne   (c version)
 * $Id$
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            01-02-1999  --  created
 * $Log$
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define SRC_VERSION 1.00
#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include <common.h>	/*  our own institute's headers */

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
 *  That's where we do the work 
 */
int main ( int argc, char **argv)
{
	int ind_read;		/* index of first non option argument */

	int verbose = FALSE;
	int pr_version = FALSE;

	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		{0,'v', NULL, 0, &verbose, TRUE, "", 
		 "Print information while calculating."}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
	char description[] = "a sample program doing nothing really";
	char nonoptarg[] = " <input_filename> [<input_filename ...]";
	extern char *program_invocation_name;		/* from errproc.h */
	extern char *program_invocation_short_name;	/* from errproc.h */
	
	errinit (argv[0]);				/* from errproc.h */
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);
	if (pr_version) { display_version(); return EXIT_SUCCESS; }
	if (ind_read == -3) return EXIT_SUCCESS;
	if (ind_read < 0 || argc - ind_read < 1) {
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

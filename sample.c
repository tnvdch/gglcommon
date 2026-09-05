/*
 *      FILE: sample
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: main program
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            16-01-1998  --  created
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define SRC_VERSION 1.00
#ifndef LIBVER
#define LIBVER ""
#endif
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
	fprintf (stderr,
		 "%s version %3.2f \n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION, LIBVER, COPYRIGHT);
	fprintf (stderr,
		 "  maintainers: E. Favey <favey@geod.ethz.ch>\n\n");
}


/*
 *  That's where we do the work 
 */
int main ( int argc, char **argv)
{
	int ind_read;		/* index of first non option argument */

	/* 
	 *  variables to hold arguments from command line parsing
	 *  and default settings for command line options 
	 */
	int verbose = FALSE;
	int pr_version = FALSE;
	int myflag = FALSE;
	int ivalue = 25;
	double dvalue = 26.4;
	char *aname = NULL;

	/* 
	 *  let's define our command line options
	 *  -h is already included in library 
	 */
	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		/* a flag type entry */
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		/* a flag type entry */
		{0,'v', NULL, 0, &verbose, TRUE, "", 
		 "Print information while calculating."},
		/* a flag type entry */
		{0,'f', NULL, 0, &myflag, TRUE, "", 
		 "My special flag that will change things."},
		/* an integer argument entry */
		{0,'i', &ivalue, CMDL_INT, NULL, 0, " <ivalue>",
		 "The ivalue you need for xyz."},
		/* an double argument entry */
		{0,'d', &dvalue, CMDL_DOUBLE, NULL, 0, " <dvalue>",
		 "A floating point argument."},
		/* an char argument entry */
		{0,'s', &aname, CMDL_CHAR, NULL, 0, " <your_name>",
		 "Please pass a string :-)"}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));

	/*
	 *  Description what the program does
	 */
	char description[] = "a sample program doing nothing really";

	/*
	 *  Information about our non option arguments
	 */
	char nonoptarg[] = " <input_filename> [<input_filename ...]";

	/* 
	 *  proper error cursing with our own command name 
	 */
	extern char *program_invocation_name;		/* from errproc.h */
	extern char *program_invocation_short_name;	/* from errproc.h */
	
	errinit (argv[0]);				/* from errproc.h */

	/*
	 * Now we process all command line options
	 */
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);
	/*
	 *  String default values are filled only down here
	 */
	if (aname == NULL) aname = strdup("weirdname");

	/*
	 *  let's evaluate if the Version flag was set
	 */
	if (pr_version) { display_version(); return EXIT_SUCCESS; }

	/*
	 *  when we got -3, the user called the -h help command
	 */
	if (ind_read == -3) return EXIT_SUCCESS;

	/*
	 *  we check if we got an error (ind_read <0) or
	 *  if we got not enough option arguments.
	 *  ind_read points to the first element of non option arguments
	 *  we check that we got at least one of them
	 */
	if (ind_read < 0 || argc - ind_read < 1) {
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
		return EXIT_FAILURE;
	}

	/*
	 *    If the verbose flag was chosen, we show the content
	 *  of our variables.
	 *    At this stage, you could call one of your subroutines
	 *  with the information and let it do the rest.
	 */
	if (verbose) {
		int i;
		fprintf (stdout, "Myflag : %d\n", myflag);
		fprintf (stdout, "ivalue : %d\n", ivalue);
		fprintf (stdout, "dvalue : %12.8f\n", dvalue);
		fprintf (stdout, "aname  : %s\n", aname);
		fprintf (stdout, "non option arguments :");
		for (i=ind_read; i<argc; i++)
			fprintf (stdout, " %s", argv[i]);
		fprintf (stdout, "\n");
	}


	/* 
	 *  clean exit upon completion of tasks
	 */
	return EXIT_SUCCESS;
}

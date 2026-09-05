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
 *            23-02-1999  --  created
 * $Log$
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <sys/stat.h>	/*  using file attributes */
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
	int every = TRUE;
	int fnum =0;

	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		{0,'n', NULL, 0, &verbose, TRUE, "", 
		 "Print number of characters read."},
		{0,'t', NULL, 0, &every, FALSE, "", 
		 "Show information only on text files."}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
	char description[] = "checks, whether a file is dos or unix text";
	char nonoptarg[] = " <input_filename> [<input_filename> ...]";
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

	for (fnum=ind_read; fnum<argc;fnum++) {
		FILE *fptxt;
		int cnum=0,c, type=1;	/* types: 0=data; 1=ascii; 2=dos */
		struct stat attrib;
		int err;

		err = stat (argv[fnum], &attrib);
		if (err!= 0) { 
			if (every) fprintf (stdout, "%s\t: %m\n", argv[fnum]); 
			continue; 
		}
		if (S_ISDIR(attrib.st_mode)) {
			if (every) fprintf (stdout, "%s\t: directory\n", argv[fnum]); 
			continue;
		}
		if (!S_ISREG(attrib.st_mode)) {
			if (every) fprintf (stdout, "%s\t: not a regular file\n", argv[fnum]); continue;
		}
		fptxt = open_text_file (argv[fnum], "r");
		if (!fptxt) continue;
		while ((c=getc(fptxt)) != EOF) {
			cnum++;
			if (isprint(c)) continue;
			if (c=='\r') {
				int second = getc(fptxt);
				cnum++;
				if (second=='\n') {
					type=2;
					break;
				}
			}
			if (c=='\n') {
				if (cnum < 2) continue;
				break;
			}
			type=0;
			break;
		}
		fclose (fptxt);
		if (cnum < 1) {
			if (every) fprintf (stdout, "%s\t: empty file\n", argv[fnum]);
		} else {
			switch (type) {
			case 1:  
				fprintf (stdout, "%s\t:", argv[fnum]);
				if (verbose) fprintf (stdout, 
						      " %3d characters read:", cnum);
				fprintf (stdout, " unix text\n"); break;
			case 2:  
				fprintf (stdout, "%s\t:", argv[fnum]);
				if (verbose) fprintf (stdout, 
						      " %3d characters read:", cnum);
				fprintf (stdout, " DOS text\n"); break;
			default: if (every) {
				fprintf (stdout, "%s\t:", argv[fnum]);
				if (verbose) fprintf (stdout, 
						      " %3d characters read:", cnum);
				fprintf (stdout, " data\n");
			}
			}
		}
	}

	return EXIT_SUCCESS;
}

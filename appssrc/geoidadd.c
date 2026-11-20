/*
 *      FILE: xyz2xyz
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: geodetic coordinate transformation
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

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>	/*  to parse the command line */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include "common.h"	/*  our own institute's headers */

/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void display_version()
{
	extern char *program_invocation_short_name;
	fprintf (stderr,
		 "%s version %s\n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION, LIBVER, COPYRIGHT);
	fprintf (stderr,
		 "  maintainers: E. Favey <favey@geod.ethz.ch>, ");
	fprintf (stderr,
		 "M. Cocard <cocard@geod.ethz.ch>\n\n");
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


static int validate_col (
	int neh_col[3],
	int neh_index[3])
{
	int i;
	for (i=0;i<3;i++) {
		neh_index[i] = -1;
	}
	for (i=0;i<3;i++) {
		if (neh_col[i] ==  4 ) neh_index[0] = i;
		if (neh_col[i] ==  5 ) neh_index[1] = i;
		if (neh_col[i] ==  6 ) neh_index[2] = i;
	}
	for (i=0;i<3;i++) {
#ifdef LIB_DEBUG
		fprintf (stderr, "neh[%d] %d ",
			 i, neh_index[i]);
#endif
		if (neh_index[i] < 0 ) {
			/* needed index not found */
			return GGLERR_INVAL;
		}
	}	
	return 0;
}


static int read_error (const char *tmp1, const char *tmp2, const unsigned long lc)
{
	extern char *program_invocation_short_name;
	if (tmp1 == tmp2)
	{
		fprintf (stderr, 
			 "%s: Unable to read value in line %ld\n",
			 program_invocation_short_name, lc);
		return GGLERR_INVAL;
	}
	return 0;
}


/*
 *  That's where we do the work 
 */
int main ( int argc, char **argv)
{
	int ind_read;
	/* for use in command line parsing */
	int verbose = FALSE;
	int pr_version = FALSE;
	char *geoidstr=NULL;
	char *inft = NULL;
	char *outft = NULL;
	int prefix=0;
	int write4=FALSE;
	int isadd=TRUE;
	/* for use in main */
	char entireline[1024];
	struct scoordcol coordstr;
	char *prefixstr;  
	char suffixstr[1024];
	int suffixsize = 1023;
	char *ptmp;
	double source[3], target[3];
	struct Grid *geoid;
	double undul;
	int g_ftyp[3], g_index[3];
	int i, err;
	unsigned long linecount;
	int in_ftyp[3], in_fcol[3];
	int out_ftyp[3], out_fcol[3];
	int in_index[3], out_index[3];

	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		/* a flag type entry */
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		{0,'v', NULL, 0, &verbose, TRUE, "", 
		 "Print information while calculating."},
		{0,'p', &prefix, CMDL_INT, NULL, 0, " <prefix_cols>",
		 "Skip <prefix_cols> columns brefore parsing values"},
		{0,'G', &geoidstr, CMDL_CHAR, NULL, 0, " <geoid>",
		 "Geoid to use (*), see $TRAFOCONF."},
		{0,'x', NULL, 0, &write4, TRUE, "", 
		 "Do not add to height, but write fourth column."},
		{0,'I', &inft, CMDL_CHAR, NULL, 0, " <input formtype>",
		 "Alternate formtype on input. Default: enh"},
		{0,'O', &outft, CMDL_CHAR, NULL, 0, " <output formtype>",
		 "Alternate formtype on output. Default: enh"}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));

	char description[] = "add or subtract geoid: h(otho)=h(ell)-undulation";
	char nonoptarg[] = "\n reads from standard in and writes to standard out\n (*) This command line switch may not be ommitted";
	
	extern char *program_invocation_name;        /* from errproc.h */
	extern char *program_invocation_short_name;  /* from errproc.h */
	
	
	errinit (argv[0]);      /* from errproc.h */
	/* default settings for command line options */
	if (strstr (program_invocation_short_name, "geoidadd") != NULL) {
		isadd = TRUE;
	} else {
		if (strstr (program_invocation_short_name, "geoidsub") != NULL) {
			isadd = FALSE;
		} else {
			fprintf (stderr,
				 "%s: invalid program name, did you rename me?\n",
				 program_invocation_short_name);
			return EXIT_FAILURE;
		}
	}
	/*
	 * Now we process all command line options
	 */
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);

	if (pr_version) { display_version(); return EXIT_SUCCESS; }
	if (ind_read == -3) return EXIT_SUCCESS;


	if (argc - optind != 0 || geoidstr == NULL) {
	  /* we do not have the necessary information or we have too much */
	  print_usage_info(FALSE, myopts, optcount, 
			   description, nonoptarg);
	  return EXIT_FAILURE;
	}
	if (inft == NULL)
		inft = strdup ("enh");
	if (outft == NULL)
		outft = strdup ("enh");

	if (verbose != FALSE) display_version();

	/*
	 *  this is where we do all the action!
	 *  read all kinds of stuff from file
	 */

	err = select_formtype(inft, in_ftyp, in_fcol);
#ifdef LIB_DEBUG
	fprintf (stderr, "%d in typ %d %d %d col %d %d %d\n", err,
		 in_ftyp[0], in_ftyp[1], in_ftyp[2],
		 in_fcol[0], in_fcol[1], in_fcol[2]);
#endif
	if (err != 0)
	{
		free(inft);
		free(outft);
		return EXIT_FAILURE;
	}
	coordstr.count = in_fcol[0] + in_fcol[1] + in_fcol[2];
	err = select_formtype(outft, out_ftyp, out_fcol);
#ifdef LIB_DEBUG
	fprintf (stderr, "%d out typ %d %d %d col %d %d %d\n", err,
		 out_ftyp[0], out_ftyp[1], out_ftyp[2],
		 out_fcol[0], out_fcol[1], out_fcol[2]);
#endif
	free(inft);
	free(outft);
    inft = NULL;
    outft = NULL;
	if (err != 0) return EXIT_FAILURE;

	geoid = select_geoid (geoidstr, g_ftyp);
	if (geoid == NULL) return EXIT_FAILURE;

	linecount = 0L;

	/* validate column types and resolve indices */
	err = validate_col (in_ftyp, in_index);
	if (err == 0)
		err = validate_col (out_ftyp, out_index);
	if (err != 0) {
		fprintf (stderr, "%s: Tell me how to resolve these columns...\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: I would guess a misconfigured config file or\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: are you trying something like latitude in meters?\n",
			 program_invocation_short_name);
		g_free (geoid);
		return EXIT_FAILURE;
	}
	err = validate_col (g_ftyp, g_index);
	if (err != 0) {
		fprintf (stderr, "%s: Tell me what style of geoid...\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: I would guess a misconfigured config file\n",
			 program_invocation_short_name);
        g_free (geoid);
		return EXIT_FAILURE;
	}

/*
 *  allocate static line space buffer
 */
	prefixstr = (char *) malloc ((prefix+2)*sizeof(char));
	coordstr.count = 3;
	coordstr.slen = 100;
	coordstr.col = (char **) malloc (coordstr.count * sizeof(char *));
	for (i=0;i<coordstr.count;i++) {
		coordstr.col[i] = (char *) malloc (coordstr.slen * sizeof(char));
	}
	while (read_line_static (stdin, entireline, 1024) != EOF) {
		linecount++;
		//if (entireline == NULL) break;
		if (entireline[0] == '#' || entireline[0] == '$' ||
		    entireline[0] == '\n'|| entireline[0] == '\r') {
			fprintf (stdout, "%s", entireline);
			continue;
		}
		err = extract_values_static (entireline, prefixstr, prefix,
					     &coordstr, suffixstr, suffixsize,
					     linecount);
		if (err != 0) {
			fprintf (stderr, "%s: invalid input line nr %ld\n",
				 program_invocation_short_name, linecount);
			fprintf (stdout, "%s", entireline);
			free(prefixstr);
			for (i=0;i<coordstr.count;i++) free (coordstr.col[i]);
			free (coordstr.col);
			g_free (geoid);
			return EXIT_FAILURE;
		}
		/* convert */
		err = 0;
		for (i=0;i<3;i++) {
			ptmp = NULL;
			source[i] = strtod(coordstr.col[in_index[i]], &ptmp);
			err = read_error (ptmp,
					  coordstr.col[in_index[i]],
					  linecount);
		}
		if (err != 0){
			free(prefixstr);
			for (i=0;i<coordstr.count;i++) free (coordstr.col[i]);
			free (coordstr.col);
			g_free (geoid);
			return EXIT_FAILURE;
		}

		/* conversion routines */
		g_inter_bilinear (geoid, source[g_index[0]], source[g_index[1]], &undul);
		target[0] = source[0];
		target[1] = source[1];
		if (write4) 
			target[2] = source[2];
		else
			target[2] = source[2] + ((isadd != FALSE ) ? undul : -undul);

		for (i=0;i<3;i++) {
			char newstr[256];
			sprintf (newstr, "%14.4f", target[i]);
			strncpy (coordstr.col[out_index[i]],newstr,
				 coordstr.slen);
		}

		if (prefix > 0) {
			fprintf (stdout, "%s", prefixstr);
		}
		for (i=0;i<3;i++) {
			fprintf (stdout, "%s ", coordstr.col[i]);
		}
		if (write4) {
			fprintf (stdout, "%14.4f ", undul);
		}
		if (strlen(suffixstr) > 0) {
			fprintf (stdout, "%s", suffixstr);
		}
		fprintf (stdout, "\n");
	}

	g_free (geoid);
	for (i=0;i<coordstr.count;i++) free (coordstr.col[i]);
	free (coordstr.col);
	free (prefixstr);
	return EXIT_SUCCESS;
}


 


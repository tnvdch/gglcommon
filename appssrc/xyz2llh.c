/*
 *      FILE: xyz2llh
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
 *            07-05-1997  --  created
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
	int xyz_col[3],
	int xyz_index[3],
	int llh_col[3],
	int llh_index[3])
{
	int i;
	for (i=0;i<3;i++) {
		xyz_index[i] = llh_index[i] = -1;
	}
	for (i=0;i<3;i++) {
		if (xyz_col[i] ==  1 ) xyz_index[0] = i;
		if (xyz_col[i] ==  2 ) xyz_index[1] = i;
		if (xyz_col[i] ==  3 ) xyz_index[2] = i;
		if (llh_col[i] ==  7 ||
		    llh_col[i] ==  9 ||
		    llh_col[i] == 11 ) llh_index[0] = i;
		if (llh_col[i] == 13 ) llh_index[0] = i;
		if (llh_col[i] ==  8 ||
		    llh_col[i] == 10 ||
		    llh_col[i] == 12 ) llh_index[1] = i;
		if (llh_col[i] == 14 ) llh_index[1] = i;
		if (llh_col[i] ==  6 ) llh_index[2] = i;
		/* h stays fixed to third col in llh, */
		/* need to add check for 13 and 14 to change */
	}
	for (i=0;i<3;i++) {
#ifdef LIB_DEBUG
		fprintf (stderr, "xyz[%d] %d llh[%d] %d ",
			 i, xyz_index[i], i, llh_index[i]);
#endif
		if (xyz_index[i] < 0 || llh_index[i] < 0) {
			/* needed index not found */
			return GGLERR_INVAL;
		}
	}	
	return 0;
}



/*
 * reincludes leading blanks 
 */ 
int with_blanks(char *s, int count)
{
  int i=count;
  while (i>0 && (s[i-1]==' ' || s[i-1]=='\n')) i--;
  return i;
}



/*
 *  That's where we do the work 
 */
int main ( int argc, char **argv)
{
	int ind_read;
	/* for use in command line parsing */
	int pr_version = FALSE;
	int verbose = FALSE;
	char *ellstr = NULL;
	char *inft = NULL;
	char *outft = NULL;
	int prefix = 0;
	int isx2e = FALSE;
	/* for use in main */
	char entireline[1024];
	char *s;  
	double source[3], target[3];
	struct ellipsoid *useell;
	int i, err;
	unsigned long linecount;
	int in_ftyp[3], in_fcol[3];
	int out_ftyp[3], out_fcol[3];
	int xyz_index[3], llh_index[3];
	
	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		/* a flag type entry */
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		/* a flag type entry */
		{0,'v', NULL, 0, &verbose, TRUE, "", 
		 "Print information while calculating."},
		{0,'e', &ellstr, CMDL_CHAR, NULL, 0, " <ellipsoid>",
		 "Use this ellipsoid for conversion on stream\n        This command line switch may not be ommitted\n        supported ellipsoids: consult the file $TRAFOCONF"},
		{0,'p', &prefix, CMDL_INT, NULL, 0, " <prefix_cols>",
		 "Skip <prefix_cols> columns brefore parsing values"},
		{0,'I', &inft, CMDL_CHAR, NULL, 0, " <input formtype>",
		 "Alternate formtype on input."},
		{0,'O', &outft, CMDL_CHAR, NULL, 0, " <output formtype>",
		 "Alternate formtype on output."}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
	
	char description[] = "ellipsoid coodinates <--> geocentric cartesian";
	char nonoptarg[] = "\n reads from standard in and writes to standard out";

	extern char *program_invocation_name;        /* from errproc.h */
	extern char *program_invocation_short_name;  /* from errproc.h */
	
	errinit (argv[0]);      /* from errproc.h */

	/* default settings for command line options */
	if (strstr (program_invocation_short_name, "xyz2llh") != NULL) {
		isx2e = TRUE;
	} else {
		if (strstr (program_invocation_short_name, "llh2xyz") != NULL) {
			isx2e = FALSE;
		} else {
			fprintf (stderr,
				 "%s: invalid program name, did you rename me?\n",
				 program_invocation_short_name);
			return EXIT_FAILURE;
		}
	}
	verbose = FALSE;
	ellstr = inft = outft = NULL;
	prefix = 0;
	/*
	 * Now we process all command line options
	 */
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);

	if (pr_version) { display_version(); return EXIT_SUCCESS; }
	if (ind_read == -3) return EXIT_SUCCESS;

	if (argc - optind != 0 || ellstr == NULL) {
		/* we do not have the necessary information or we have too much */
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
		return EXIT_FAILURE;
	}
	if (inft == NULL) {
		if (isx2e != FALSE) {
			inft = strdup ("xyz");
		} else {
			inft = strdup ("llh");
		}
	}
	if (outft == NULL) {
		if (isx2e != FALSE) {
			outft = strdup ("llh");
		} else {
			outft = strdup ("xyz");
		}
	}

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
	    inft = NULL;
	    free(outft);
	    outft = NULL;
		return EXIT_FAILURE;
	}

	err = select_formtype(outft, out_ftyp, out_fcol);
#ifdef LIB_DEBUG
	fprintf (stderr, "%d out typ %d %d %d col %d %d %d\n", err,
		 out_ftyp[0], out_ftyp[1], out_ftyp[2],
		 out_fcol[0], out_fcol[1], out_fcol[2]);
#endif
	if (err != 0)
	{
		free(inft);
	    inft = NULL;
	    free(outft);
	    outft = NULL;
		return EXIT_FAILURE;
	}
	useell = select_ellipsoid (ellstr);
	if (useell == NULL)
	{
		free(inft);
	    inft = NULL;
	    free(outft);
	    outft = NULL;
		return EXIT_FAILURE;
	}
	free(inft);
	inft = NULL;
    free(outft);
	outft = NULL;

	linecount = 0L;

	/* validate column types and resolve indices */
	if (isx2e != FALSE)
		err = validate_col (in_ftyp, xyz_index, out_ftyp, llh_index);
	else
		err = validate_col (out_ftyp, xyz_index, in_ftyp, llh_index);
	if (err != 0) {
		fprintf (stderr, "%s: Tell me how to resolve these columns...\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: I would guess a misconfigured config file or\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: are you trying something like latitude in meters?\n",
			 program_invocation_short_name);

		return EXIT_FAILURE;
	}

	while (read_line_static (stdin, entireline, 1024) != EOF) {
		int count;
		linecount++;
		//if (entireline == NULL) break;
		if (entireline[0] == '#' || entireline[0] == '$' ||
		    entireline[0] == '\n'|| entireline[0] == '\r') {
			fprintf (stdout, "%s", entireline);
			continue;
		}

		if (prefix >= strlen(entireline)) { /* error in inputline */
			fprintf (stderr, "%s: invalid input in line %ld (prefix longer than line)\n",
				 program_invocation_short_name, linecount);
			return EXIT_FAILURE;
		}

		s=entireline+prefix; /* position to start of time info in string */

		/* convert */

		err = 0;
		if (isx2e != FALSE) {
			struct Angle lat, lon;
			char newstr[256];
			if (sscanf(s," %lf %lf %lf %n", &(source[xyz_index[0]]),
				   &(source[xyz_index[1]]), &(source[xyz_index[2]]),
				   &count) != 3) {
				fprintf (stderr, 
					 "%s:error in converting line nr %ld\n",
					 program_invocation_short_name, linecount);
				return EXIT_FAILURE;
			}
			xyz_to_llh ((*useell), source,
				    &target[0], &target[1], &target[2]);
			entireline[prefix]='\0';  /* clip prefix in entire line */
			fprintf(stdout,"%s",entireline); /* write prefix */
			lat = afromrad (target[0]);
			lon = afromrad (target[1]);
			for (i=0;i<3;i++) {
				switch (out_ftyp[i]) {
				case  7: sprintf (newstr, "%14.10f", lat.deg); break;
				case  8: sprintf (newstr, "%14.10f", lon.deg); break;
				case  9: sprintf (newstr, "%14.12f",lat.rad); break;
				case 10: sprintf (newstr, "%14.12f",lon.rad); break;
				case 11: sprintf (newstr, "%14.10f",lat.grad); break;
				case 12: sprintf (newstr, "%14.10f",lon.grad); break;
				case 13: sprintf (newstr, "%s", admsstr(lat,0)); break;
				case 14: sprintf (newstr, "%s", admsstr(lon,0)); break;
				case  6: sprintf (newstr, "%14.4f", target[2]); break;
#ifdef WITH_DEBUG
				default: fprintf (stderr, "%s: internal error\n",
						  program_invocation_short_name);
					 abort();
#endif
				}
				fprintf (stdout, " %s", newstr);
			}
		} else {
			struct Angle lat, lon;

			if ((in_ftyp[0]< 13) && (in_ftyp[1] < 13)) {
				if (sscanf(s," %lf %lf %lf %n", 
					   &(source[llh_index[0]]),
					   &(source[llh_index[1]]), 
					   &(source[llh_index[2]]),
					   &count) != 3) {
					fprintf (stderr, 
						 "%s:error in converting line nr %ld\n",
						 program_invocation_short_name, linecount);
					return EXIT_FAILURE;
				}
		
				for (i=0;i<3;i++) {
					switch (in_ftyp[llh_index[i]]) {
					case  7: lat = afromdeg (source[i]); 
						source[i] = lat.rad; break;
					case  8: lon = afromdeg (source[i]);
						source[i] = lon.rad; break;
					case  9: lat = afromrad (source[i]);
						source[i] = lat.rad; break;
					case 10: lon = afromrad (source[i]);
						source[i] = lon.rad; break;
					case 11: lat = afromgrad(source[i]);
						source[i] = lat.rad; break;
					case 12: lon = afromgrad(source[i]);
						source[i] = lon.rad; break;
					case  6: break;
#ifdef WITH_DEBUG
					default: fprintf (stderr, "%s: internal error\n",
							  program_invocation_short_name);
					abort();
#endif
					}
				}
			} else {
				/* the dms case */
				int deg[2], min[2];
				double sec[2];
				if (sscanf(s," %d %d %lf %d %d %lf %lf %n",
					   &(deg[llh_index[0]]), 
					   &(min[llh_index[0]]), 
					   &(sec[llh_index[0]]),
					   &(deg[llh_index[1]]), 
					   &(min[llh_index[1]]), 
					   &(sec[llh_index[1]]), 
					   /* height in the dms case */
					   &(source[llh_index[2]]),
					   &count) != 3) {
					fprintf (stderr, 
						 "%s:error in converting line nr %ld\n",
						 program_invocation_short_name, linecount);
					return EXIT_FAILURE;
				}
				for (i=0;i<2;i++) if (in_ftyp[i] == 13) {
					lat = afromdmsverb (deg[i], min[i], sec[i]);
					source[i] = lat.rad;
				} else {
					lon = afromdmsverb (deg[i], min[i], sec[i]);
					source[i] = lon.rad;
				} /* for */
			}
#ifdef LIB_DEBUG
			fprintf (stderr, "%14.9f %14.9f %14.9f\n", source[0], source[1], source[2]);
#endif
			llh_to_xyz ((*useell), source[0], source[1], source[2], target);
			entireline[prefix]='\0';  /* clip prefix in entire line */
			fprintf(stdout,"%s",entireline); /* write prefix */
			for (i=0;i<3;i++) {
				fprintf (stdout, " %14.4f", target[xyz_index[i]]);
			}
		}
		s+=with_blanks(s,count);     /* position to start of postfix (always contains at least \n */
		fprintf(stdout,"%s",s);          /* write postfix */

	} /* while <STDIN> */

	return EXIT_SUCCESS;
}


 


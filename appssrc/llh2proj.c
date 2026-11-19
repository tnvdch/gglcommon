/*
 *      FILE: llh2proj
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

/*
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void display_usage (const int verbose, const int isp2e)
{
	extern char *program_invocation_short_name;
	fprintf (stderr,
		 "usage:\n   %s [-Vhv] -P <projection> [-p <prefix_cols>]\n" ,
		 program_invocation_short_name);
	fprintf (stderr,
		 "           [-I <input formtype>] [-O <output formtype>]\n\n");
	fprintf (stderr,"  reads from standard in and writes to standard out\n\n");
	
	if (verbose != FALSE)
	{
		char *indef, *outdef;
		indef = outdef = NULL;
		
		if (isp2e != FALSE) {
			indef = strdup ("enh"); outdef = strdup ("llh");
		} else {
			indef = strdup ("llh"); outdef = strdup ("enh");
		}

		fprintf (stderr,"  -V    Version string\n");
		fprintf (stderr,"  -h    This help message\n");
		fprintf (stderr,"  -v    verbose mode, log conversion progress to stderr\n");
		fprintf (stderr,"  -P    Use this projection for conversion on stream\n");
		fprintf (stderr,"        This command line switch may not be ommitted\n");
		fprintf (stderr,"        supported projections: consult the file $TRAFOCONF\n");
		fprintf (stderr,"  -p    Skip <prefix_cols> columns brefore parsing values\n");
		fprintf (stderr,"  -I    Alternate formtype on input. Default: %s\n", indef);
		fprintf (stderr,"  -O    Alternate formtype on output. Default: %s\n\n", outdef);
		fprintf (stderr,"        supported formtypes: consult the file $TRAFOCONF\n");
		free (indef); free (outdef);
	}
	else
	{
		fprintf (stderr,
			 "  -h    Prints more detailed help on each option\n\n");      
	}
}


static int validate_col (
	int neh_col[3],
	int neh_index[3],
	int llh_col[3],
	int llh_index[3])
{
	int i;
	for (i=0;i<3;i++) {
		neh_index[i] = llh_index[i] = -1;
	}
	for (i=0;i<3;i++) {
		if (neh_col[i] ==  4 ) neh_index[0] = i;
		if (neh_col[i] ==  5 ) neh_index[1] = i;
		if (neh_col[i] ==  6 ) neh_index[2] = i;
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
#ifdef WITH_DEBUG
		fprintf (stderr, "neh[%d] %d llh[%d] %d ",
			 i, neh_index[i], i, llh_index[i]);
#endif
		if (neh_index[i] < 0 || llh_index[i] < 0) {
			/* needed index not found */
			return GGLERR_INVAL;
		}
	}	
	return 0;
}


static int read_error (const char *tmp1, const char *tmp2, const unsigned long lc)
{
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
	int arg_char;		/* argument character to process */
	/* for use in command line parsing */
	int verbose;
	char *projstr;
	char *inft;
	char *outft;
	int prefix;
	int isp2e;
	/* for use in main */
	char *entireline;
	struct coordcol coordstr;
	char *prefixstr;  
	char *suffixstr;
	char *ptmp;
	double source[3], target[3];
	int i, err;
	unsigned long linecount;
	int in_ftyp[3], in_fcol[3];
	int out_ftyp[3], out_fcol[3];
	int neh_index[3], llh_index[3];
	
	extern char *program_invocation_name;        /* from errproc.h */
	extern char *program_invocation_short_name;  /* from errproc.h */
	
	extern int opterr;      /* from unistd.h: let's do our own cursing */
	
	errinit (argv[0]);      /* from errproc.h */
	opterr = 0;
	/* default settings for command line options */
	if (strstr (program_invocation_short_name, "proj2llh") != NULL) {
		isp2e = TRUE;
	} else {
		if (strstr (program_invocation_short_name, "llh2proj") != NULL) {
			isp2e = FALSE;
		} else {
			fprintf (stderr,
				 "%s: invalid program name, did you rename me?\n",
				 program_invocation_short_name);
			return EXIT_FAILURE;
		}
	}
	verbose = FALSE;
	projstr = inft = outft = NULL;
	prefix = 0;
	/*
	 * Now we process all command line options
	 */
	while ((arg_char = getopt (argc, argv, "VhvP:p:O:I:")) != -1)
		switch (arg_char)
		{
		case 'h': display_usage(TRUE, isp2e);	/* -h   the help option */
			/* fallthrough */
		case 'V': display_version();		/* -v   the version string */
			return EXIT_SUCCESS;
		case 'v': verbose = TRUE;    	        break;
		case 'P': projstr = strdup (optarg);	break;
		case 'p':
		{
			char *ptmp = NULL;
			prefix = strtol(optarg, &ptmp, 10);
			if (ptmp == optarg || prefix < 0)
			{
				fprintf (stderr, 
					 "%s: Option `-l' requires positive integer argument.\n",
					 program_invocation_short_name);
				display_usage(FALSE, isp2e);
				return EXIT_FAILURE;
			}
		}
		break;
		case 'I': inft = strdup (optarg);  break;
		case 'O': outft = strdup (optarg); break;
		
		case '?':                              /* oops, an unknown option */
			if (isprint (optopt))
				fprintf (stderr, "%s: Unknown option `-%c'.\n",
					 program_invocation_short_name, optopt);
			else
				fprintf (stderr,
					 "%s: Unknown option character `\\x%x'.\n",
					 program_invocation_short_name, optopt);
			display_usage(FALSE, isp2e);
			return EXIT_FAILURE;
			
		default:              /* actually, this should never happen */
			abort ();
		}

	if (argc - optind != 0 || projstr == NULL) {
		/* we do not have the necessary information or we have too much */
		display_usage(FALSE, isp2e);
		return EXIT_FAILURE;
	}
	if (inft == NULL) {
		if (isp2e != FALSE) {
			inft = strdup ("enh");
		} else {
			inft = strdup ("llh");
		}
	}
	if (outft == NULL) {
		if (isp2e != FALSE) {
			outft = strdup ("llh");
		} else {
			outft = strdup ("enh");
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
	if (err != 0)
	{
		free(inft);
		free(outft);
		return EXIT_FAILURE;
	}

	free(inft);
	free(outft);
	inft = NULL;
	outft = NULL;

	/* here we would have to switch (projection) */
	/* verify the projection type, at the time being only swiss */
	if (strstr(projstr, "swiss") == NULL)
	{
		return EXIT_FAILURE;
	}

	linecount = 0L;

	/* validate column types and resolve indices */
	if (isp2e != FALSE)
		err = validate_col (in_ftyp, neh_index, out_ftyp, llh_index);
	else
		err = validate_col (out_ftyp, neh_index, in_ftyp, llh_index);
	if (err != 0) {
		fprintf (stderr, "%s: Tell me how to resolve these columns...\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: I would guess a misconfigured config file or\n",
			 program_invocation_short_name);
		fprintf (stderr, "%s: are you trying something like latitude in meters?\n",
			 program_invocation_short_name);
		return EXIT_FAILURE;
	}

	while (read_line (stdin, &entireline) != EOF) {
		linecount++;
		if (entireline == NULL) break;
		if (entireline[0] == '#' || entireline[0] == '$' ||
		    entireline[0] == '\n'|| entireline[0] == '\r') {
			fprintf (stdout, "%s", entireline);
			free (entireline);
			continue;
		}
		err = extract_values (entireline, &prefixstr, prefix,
				      &coordstr, &suffixstr, linecount);
		if (err != 0) {
			fprintf (stderr, "%s: invalid input line nr %ld\n",
				 program_invocation_short_name, linecount);
			fprintf (stdout, "%s", entireline);
			free (entireline);
			return EXIT_FAILURE;
		}
		free (entireline);
		/* convert */
		err = 0;
		if (isp2e != FALSE) {
			struct Angle lat, lon;
			char newstr[256];
			for (i=0;i<3;i++) {
				ptmp = NULL;
				source[i] = strtod(coordstr.col[neh_index[i]], &ptmp);
				err = read_error (ptmp,
						  coordstr.col[neh_index[i]],
						  linecount);
				free (coordstr.col[neh_index[i]]);
			}
			if (err != 0) return EXIT_FAILURE;
			/* here we would have to switch (projection) */
			ch_miltoell (source[0], source[1],
				    &target[0], &target[1]);
			target[2] = source[2];
			lat = afromrad (target[0]);
			lon = afromrad (target[1]);
			for (i=0;i<3;i++) {
				switch (out_ftyp[i]) {
				case  7: sprintf (newstr, "%14.9f", lat.deg); break;
				case  8: sprintf (newstr, "%14.9f", lon.deg); break;
				case  9: sprintf (newstr, "%14.12f",lat.rad); break;
				case 10: sprintf (newstr, "%14.12f",lon.rad); break;
				case 11: sprintf (newstr, "%14.9f",lat.grad); break;
				case 12: sprintf (newstr, "%14.9f",lon.grad); break;
				case 13: sprintf (newstr, "%s", admsstr(lat,0)); break;
				case 14: sprintf (newstr, "%s", admsstr(lon,0)); break;
				case  6: sprintf (newstr, "%14.3f", target[2]); break;
#ifdef WITH_DEBUG
				default: fprintf (stderr, "%s: internal error\n",
						  program_invocation_short_name);
					 abort();
#endif
				}
				coordstr.col[llh_index[i]] = strdup (newstr);
			}
		} else {
			struct Angle lat, lon;
			char newstr[256];

			if ((in_ftyp[0]< 13) && (in_ftyp[1] < 13)) {
				for (i=0;i<3;i++) {
					ptmp = NULL;
					source[i] = strtod(coordstr.col[llh_index[i]], &ptmp);
					err = read_error (ptmp,
							  coordstr.col[llh_index[i]],
							  linecount);
					free (coordstr.col[llh_index[i]]);
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
				int deg, min;
				double sec;
				err = 0;
				for (i=0;i<2;i++) {
					ptmp = NULL;
					deg = strtol(coordstr.col[llh_index[i]*3], &ptmp, 10);
					err = read_error (ptmp,
							  coordstr.col[llh_index[i]*3],
							  linecount);
					ptmp = NULL;
					min = strtol(coordstr.col[llh_index[i]*3+1], &ptmp, 10);
					err = read_error (ptmp,
							  coordstr.col[llh_index[i]*3+1],
							  linecount);
					ptmp = NULL;
					sec = strtod(coordstr.col[llh_index[i]*3+2], &ptmp);
					err = read_error (ptmp,
							  coordstr.col[llh_index[i]*3+2],
							  linecount);
					if (in_ftyp[i] == 13) {
						lat = afromdmsverb (deg, min, sec);
						source[i] = lat.rad;
					} else {
						lon = afromdmsverb (deg, min, sec);
						source[i] = lon.rad;
					}
					free (coordstr.col[llh_index[i]*3]);
					free (coordstr.col[llh_index[i]*3+1]);
					free (coordstr.col[llh_index[i]*3+2]);
				}
				/* height in the dms case */
				ptmp = NULL;
				source[2] = strtod(coordstr.col[llh_index[2]*3], &ptmp);
				err = read_error (ptmp, coordstr.col[llh_index[2]*3],
						  linecount);
				free (coordstr.col[llh_index[2]*3]);
			}
			if (err != 0) return EXIT_FAILURE;
#ifdef LIB_DEBUG
			fprintf (stderr, "%14.9f %14.9f %14.9f\n", source[0], source[1], source[2]);
#endif
			/* here we would have to switch (projection) */
			ch_elltomil (source[0], source[1], &target[0], &target[1]);
			target[2] = source[2];
			for (i=0;i<3;i++) {
				sprintf (newstr, "%14.4f", target[i]);
				coordstr.col[neh_index[i]] = strdup (newstr);
			}
		}
		if (prefixstr != NULL) {
			fprintf (stdout, "%s", prefixstr);
			free (prefixstr);
		}
		for (i=0;i<3;i++) {
			fprintf (stdout, "%s ", coordstr.col[i]);
			free (coordstr.col[i]);
		}
		if (suffixstr != NULL) {
			fprintf (stdout, "%s", suffixstr);
			free (suffixstr);
		}
		fprintf (stdout, "\n");
	}

	return EXIT_SUCCESS;
}


 


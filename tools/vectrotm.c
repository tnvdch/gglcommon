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

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include <common.h>	/*  our own institute's headers */

#define RHO (180.0/M_PI)
/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
static void display_version()
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

/***********************************************************************/

static void enu2nel (long double vect[3])
{
	long double b = vect[0];
	vect[0] = vect[1];
	vect[1] = b;
	vect[2] = -vect[2];
}

/***********************************************************************/

static void calcivect (
	const long double v1[3], const long double v2[3], long double e[3][3])
{
	long double t1[3];
	lv3_norm(v1, e[0]);         lv3_norm(v2, t1); 
	lv3_cross(e[0], t1, e[2]);  lv3_cross(e[2], e[0], e[1]);
	lv3_norm(e[1], e[1]);         lv3_norm(e[2], e[2]); 
}

/***********************************************************************/

static void calctmat (		/* rotation matrix system 1 -> system 2 */
	const long double e1[3][3],
	const long double e2[3][3],
	long double tmat[3][3])
{
	int i,j;
	for (i=0;i<3;i++) for (j=0;j<3;j++)
		tmat[i][j] = lv3_scalar(e1[i], e2[j]);
}

/***********************************************************************/

static void detangle(
	const long double rotmat[3][3], double *roll, double *pitch, double *az)
{
	*pitch = asin(-rotmat[2][0]);
	*az = atan2 (rotmat[1][0], rotmat[0][0]);
	*roll = atan2 (rotmat[2][1], rotmat[2][2]);
}

/***********************************************************************/

#ifdef KUNST
static void dorot (
	long double src[3], long double trg[3], double roll, double pitch, double yaw)
{
	long double tmp[3];
	lrotate (src, 1, -roll/RHO, trg);
 	lrotate (trg, 2, -pitch/RHO, tmp);
 	lrotate (tmp, 3, -yaw/RHO, trg);
}
#endif

/***********************************************************************/


static void vect2angle (
	long double v1[3], long double v2[3], 
	long double vc1[3], long double vc2[3],
	double *roll, double *pitch, double *yaw)
{
	const long double unit[3][3] = {{1.0, 0.0, 0.0},
					{0.0, 1.0, 0.0},
					{0.0, 0.0, 1.0}};
	double r, p, y;
	long double ea[3][3], ec[3][3], rotmat[3][3], totmat[3][3];
	long double calmat[3][3], calmati[3][3];
/* 	long double a, rv[3]; */

	calcivect(vc1, vc2, ec);
	calcivect(v1, v2, ea);
	calctmat(unit, ec, calmat);
	lm3_trans (calmat, calmati);
	calctmat(unit, ea, rotmat);
	lm3m3_mult(rotmat, calmati, totmat);

	detangle(totmat, &r, &p, &y);
	if (y<0.0) y+= 2.0*M_PI;
	*yaw = y;
	*pitch = p;
	*roll = r;
}

/***********************************************************************/

int vectrotm (
	FILE *fpin, FILE *fpout, FILE *fperr,
	const int doppel, const int verbose)
{
	/* calibration vectors defining the aircraft system */
   	long double vc1[3] = { 0.296711,  3.921529,  0.127902};
   	long double vc2[3] = {-0.233085, -5.351753,  2.820545};
   	long double vd1[3] = {-2.713501,  3.462421,  0.363382};
   	long double vd2[3] = { 2.814676,  3.462421,  0.363382};

	long double zeit=0.0, frnt[3], tail[3], left[3], rght[3];
	double defr, defp, defy,defrd,defpd,defyd,yawdiff;
#ifndef KUNST
	int err;
#endif

	enu2nel (vc1); enu2nel (vc2);
	enu2nel (vd1); enu2nel (vd2);

#ifndef KUNST
	overread_comment (fpin, NULL, '#');
	while((err=fscanf(fpin," %Lf %Lf %Lf %Lf %Lf %Lf %Lf", &zeit,
			  &(frnt[0]), &(frnt[1]), &(frnt[2]),
			  &(tail[0]), &(tail[1]), &(tail[2])) ) != EOF) {
		enu2nel (frnt); enu2nel (tail);
		if (doppel) {
			err=fscanf(fpin," %Lf %Lf %Lf %Lf %Lf %Lf",
				   &(left[0]), &(left[1]), &(left[2]),
				   &(rght[0]), &(rght[1]), &(rght[2]));
			enu2nel (left); enu2nel (rght);
		}
#endif
#ifdef KUNST
 	zeit=0.0;
 	while (zeit < 360.0) {
 		dorot (vc1, frnt, -30.0, 3.0, zeit);
 		dorot (vc2, tail, -30.0, 3.0, zeit);
		if (doppel) {
			dorot (vd1, left, -30.0, 3.0, zeit);
			dorot (vd2, rght, -30.0, 3.0, zeit);
		}
#endif		
/* from frnt and tail */
		vect2angle (frnt, tail, vc1, vc2, &defr, &defp, &defy);
		if (doppel) fprintf (fpout, "W12 ");
		fprintf (fpout, "%12.8Lf %12.8f %12.8f %12.8f\n", zeit,
			 defr*RHO, defp*RHO, defy*RHO);

		if (doppel) {
/* from left and rght */
			vect2angle (left, rght, vd1, vd2, &defrd, &defpd, &defyd);
			fprintf (fpout, "W34 %12.8Lf %12.8f %12.8f %12.8f\n", zeit,
				 defrd*RHO, defpd*RHO, defyd*RHO);
/* difference */
			yawdiff = defy-defyd;
			if (yawdiff < -M_PI) yawdiff += 2*M_PI;
			if (yawdiff >  M_PI) yawdiff -= 2*M_PI;
			fprintf (fpout, "diff %12.8Lf %12.8f %12.8f %12.8f\n", zeit,
				 (defr-defrd)*RHO, (defp-defpd)*RHO, yawdiff*RHO);
		}
#ifdef KUNST
 		zeit += 1.0;
#endif
	}
	return 0;
}

/***********************************************************************/


int main ( int argc, char **argv)
{
	int ind_read;		/* index of first non option argument */

	int verbose = FALSE;
	int doppel = TRUE;
	int pr_version = FALSE;
	FILE *fpvec;

	struct mainoptions myopts[] =			/* from cmdline.h */
	{
		{0,'V', NULL, 0, &pr_version, TRUE, "",
		 "Version string."},
		{0,'v', NULL, 0, &verbose, TRUE, "", 
		 "Print information while calculating."},
		{0,'2', NULL, 0, &doppel, FALSE, "", 
		 "input only 2 vectors for orientation."}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
	char description[] = "compute rotation matrix and angles from vectors";
	char nonoptarg[] = " [<vector_filename>]";
	extern char *program_invocation_name;		/* from errproc.h */
	extern char *program_invocation_short_name;	/* from errproc.h */
	
	errinit (argv[0]);				/* from errproc.h */
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);
	if (pr_version) { display_version(); return EXIT_SUCCESS; }
	if (ind_read == -3) return EXIT_SUCCESS;
	if (ind_read < 0 || argc - ind_read < 0) {
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
		return EXIT_FAILURE;
	}
	if (argc - ind_read == 1) {
		fpvec = open_text_file (argv[ind_read], "r");
		if (!fpvec) return EXIT_FAILURE;
	} else {
		if (isatty(STDIN_FILENO)) {
			fprintf (stderr, 
				 "You don't wanna type in the data, do you?\n\n");
			print_usage_info(FALSE, myopts, optcount, 
					 description, nonoptarg);
			return EXIT_FAILURE;
		}
		fpvec = stdin;
	}

	vectrotm (fpvec, stdout, stderr, doppel, verbose);

	return EXIT_SUCCESS;
}

#undef RHO

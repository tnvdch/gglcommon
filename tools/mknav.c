/*
 *      FILE: $Source$
 *
 *   PROJECT: GGL GPSkin Applications
 *
 *   PURPOSE: create nav file for stuttgart
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Etienne Favey
 * $Id$
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            18-06-1999  --  created, linear approximation
 * $Log$
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1999 Geodesy and Geodynamics Lab, ETH Zurich"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <common.h>

void display_version() 
{
	extern char *program_invocation_short_name;
	extern char *common_version;
	fprintf (stderr, "%s version %s\n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION, 
		 common_version, COPYRIGHT);
	fprintf (stderr, "  maintainers: E. Favey <favey@geod.ethz.ch>\n\n");
}

static void print_version(FILE *fpout, char **argv, int argc)
{
	extern char *program_invocation_short_name;
	extern char *common_version;
	int i;
	fprintf (fpout, "# %s version %s\n# %s\n# %s\n#",
		 program_invocation_short_name, SRC_VERSION, 
		 common_version, COPYRIGHT);
	for (i=0;i<argc;i++) {
	  fprintf (fpout, " %s", argv[i]);
	}
	fprintf (fpout, "\n");
}

/***************************************************************************/
#define GETMAX(a,b) ((a>b)?a:b)
static double detmax_startt (
	struct data_window *gps,
	struct data_window *roll,
	struct data_window *pitch,
	struct data_window *yaw)
{
	double result = -1.0;
	if (gps->m != NULL) 
		result = GETMAX(result,gps->m->mat[0][0]);
	if (roll->m != NULL) 
		result = GETMAX(result,roll->m->mat[0][0]);
	if (pitch->m != NULL) 
		result = GETMAX(result,pitch->m->mat[0][0]);
	if (yaw->m != NULL) 
		result = GETMAX(result,yaw->m->mat[0][0]);
	return result;
}

/***************************************************************************/

int read_coord(FILE *fp, struct matrix *m, int row)
{
	/* reads time type 0, X, Y, Z */
	struct ttime t;
	char s[1024];
	int err, gpsweek;
	double gpssec;

	if (read_line_static(fp, s, 1024)==EOF) return EOF;
	err = sscanf (s, "%d %d %d %d %d %lf %lf %lf %lf",
		      &t.year, &t.month, &t.day, &t.hour, &t.min, &t.sec,
		      &(m->mat[row][1]), &(m->mat[row][2]), &(m->mat[row][3]));

	if (err < 9) return 0;

	/* in this program, we calculate in secinweek */
	datmjd (&t);
	mjd_to_gps(t.mjd, &gpssec, &gpsweek);
	m->mat[row][0] = gpssec;
#ifdef APP_DEBUG
	fprintf (stderr, "gps_time: %12.8f sec\n", gpssec);
#endif

	return 1;
}

/***************************************************************************/

int read_imufile(FILE *fp, double imu[5])
{
	/* reads GPSweek, secinweek, roll, pitch, yaw */
	char s[1024];
	int err;

	if (overread_comment (fp, NULL, '#') == EOF) return EOF;
	if (read_line_static(fp, s, 1024)==EOF) return EOF;
	err = sscanf(s, "%lf %lf %lf %lf %lf",
		     &imu[0], &imu[1], &imu[2], &imu[3], &imu[4]);
	if (err!=5) return 0;
	return 1;
}

/***************************************************************************/

int read_drift(FILE *fp, struct matrix *m, int row)
{
	/* reads GPSweek, secinweek, offsetvalue */
	char s[1024];
	int err;

	if (read_line_static(fp, s, 1024)==EOF) return EOF;
	err = sscanf(s, "%lf %lf %lf",
		     &m->mat[row][1], &m->mat[row][0], &m->mat[row][2]);
	if (err!=3) return 0;
	return 1;
}

/***************************************************************************/

static void det_laser_pos (
	const double xyzantenna[3], /* GPS geocentric coordinates */
	const double antairframe[3],  /* GPS in aircraft frame */
	const double lasairframe[3],  /* laser in aircraft frame */
	const double roll,  /* platform [radians] */
	const double pitch, /* platform [radians] */
	const double yaw,   /* platform [radians] */
	const struct ellipsoid ell,
	double xyzlaser[3])
{
	double lat, lon, hgt;
	double antlas_af[3]; /* antenna to laser in aircraft frame */
	double antlas_tf[3]; /* antenna to laser in topo frame */
	double antlas_gf[3]; /* antenna to laser in geocentric frame */
	double antlas_tmp[3];

	/* extract latitude and longitude for later calculation */
	xyz_to_llh (ell, xyzantenna, &lat, &lon, &hgt);

	/* calculate offset vector in air frame */
	v3_sub(lasairframe, antairframe, antlas_af);

	/* rotate offset vector to topo frame */
	rotate (antlas_af,  1, -roll,  antlas_tmp);
	rotate (antlas_tmp, 2, -pitch, antlas_tmp);
	rotate (antlas_tmp, 3, -yaw,   antlas_tf);

	/* rotate offset vector to earth frame */
	rotate (antlas_tf,  2, (lat + M_PI/2), antlas_tmp);
	rotate (antlas_tmp, 3, -lon, antlas_gf);

	/* apply offset to GPS position to get laser position */
	v3_add(xyzantenna, antlas_gf, xyzlaser);
}

/***************************************************************************/

static int gps_interpol (struct data_window *coord, const double t, double pos[3])
{
#define GCM coord->m->mat
	extern char *program_invocation_short_name;
	int i,j;

	if ((i=time_index_data_window(coord, t)) == -1) {
		window_err (-4, coord->fname); return -1; 
	}
	for (j=0;j<3;j++) {
		if (d1ipol_linear (GCM[i][0], GCM[i+1][0], GCM[i][j+1], GCM[i+1][j+1],
			       t, &(pos[j]))!=0) {
			fprintf (stderr, 
				 "%s: interpolation error, division by zero\n",
				 program_invocation_short_name);
			return -1;
		}
	}
	return 0;
#undef GCM
}

/***************************************************************************/

static int angle_interpol (
	struct data_window *drift, const double t, double *result)
{
#define GCM drift->m->mat
	extern char *program_invocation_short_name;
	int i;

	if ((i=time_index_data_window(drift, t)) == -1) {
		window_err (-4, drift->fname); return -1; 
	}
	if (d1ipol_linear (GCM[i][0], GCM[i+1][0], GCM[i][2], GCM[i+1][2],
			   t, result)!=0) {
		fprintf (stderr, 
			 "%s: interpolation error, division by zero\n",
			 program_invocation_short_name);
		return -1;
	}
#ifdef APP_DEBUG
	fprintf (stderr, "%s: %12.6f %12.6f %12.6f %12.6f = %12.6f %12.6f\n",
		 drift->fname,GCM[i][0],GCM[i+1][0],GCM[i][2],GCM[i+1][2],
		 t,*result);
#endif
	return 0;
#undef GCM
}

/***************************************************************************/

static void angle_drift_corr (
	double imu[5], double rollcorr, double pitchcorr, double yawcorr)
{
	imu[2] += rollcorr;
	imu[3] += pitchcorr;
	imu[4] += yawcorr;
}

/***************************************************************************/
/*
static void yaw_drift_corr (double imu[5])
{
 	double polycoeff[3] = {0.251084, -6.38164, 36.463470+1.86402329};
/ * 	double polycoeff[2] = {-0.891406, 6.571415}; * /
	double hour, day, corr;

	hour = modf( gps_to_mjd(imu[1], (int) imu[0]), &day) * 24.0;
	corr = polycoeff[2] + (polycoeff[1] + polycoeff[0]*hour)*hour;
#ifdef APP_DEBUG
	fprintf (stderr, "%10.7f %12.8f\n", hour, corr);
#endif
	imu[4] += corr;
}
*/
/***************************************************************************/
/*
static void linear_drift_corr (double imu[5])
{
#define TOFFSET (35400.0)
	double rollcoeff[2] = { 5.996282e-05,  1.298696};
	double pitchcoeff[2]= { 1.805594e-05, -0.7019722};
	double yawcoeff[2]  = {-0.0001478801, -1.794834};
	double secofday, day, corr;

	secofday = modf( gps_to_mjd(imu[1], (int) imu[0]), &day)*86400.0 - TOFFSET;
	corr = rollcoeff[1] + rollcoeff[0]*secofday;
	imu[2] += corr;

	corr = pitchcoeff[1] + pitchcoeff[0]*secofday;
	imu[3] += corr;

	corr = yawcoeff[1] + yawcoeff[0]*secofday;
	imu[4] += corr;

#undef TOFFSET
}
*/
/***************************************************************************/
#define RHO (180.0/M_PI)

static void mknav1 (
	const char *gpsfname, 
	const char *rollfn, const char *pitchfn, const char *yawfn,
	FILE *fpimu, FILE *fpout, const int dooff,
	const double halfwin, const int quiet)
{
	extern char *program_invocation_short_name;
	struct data_window coord, rd, pd, yd;
	int err=0;
	double imu[5];  /* 0=GPSweek, 1=secofweek, 2=roll, 3=pitch, 4=yaw */
	double gpspos[3], laspos[3];
	double antairframe[3] = {-5.4097, 0.0574, -1.5296};
	double lasairframe[3] = {-0.0455, 0.0003, -0.1300};
	struct ellipsoid wgs84 = {6378137.000, 0.00669438};
	double startt;

	/* if this stays NULL, this file was not supplied */
	coord.m = NULL; rd.m=NULL; pd.m=NULL; yd.m=NULL;

	if (quiet == FALSE && dooff != FALSE) {
	  fprintf (fpout, "# using antenna position: %7.4f %7.4f %7.4f\n",
		   antairframe[0], antairframe[1], antairframe[2]);
	  fprintf (fpout, "# using  laser  position: %7.4f %7.4f %7.4f\n",
		   lasairframe[0], lasairframe[1], lasairframe[2]);
	}

	if (gpsfname != NULL)
		if ((err = init_data_window(halfwin, 20, 4, gpsfname, 
					    &coord, read_coord)) != 0) {
			window_err(err, gpsfname); exit(err);
		} 
	if (rollfn != NULL)
		if ((err = init_data_window(halfwin, 20, 3, rollfn, 
					    &rd, read_drift)) != 0) {
			window_err(err, rollfn); exit(err);
		}
	if (pitchfn != NULL)
		if ((err = init_data_window(halfwin, 20, 3, pitchfn, 
					    &pd, read_drift)) != 0) {
			window_err(err, pitchfn); exit(err);
		}
	if (yawfn != NULL)
		if ((err = init_data_window(halfwin, 20, 3, yawfn, 
					    &yd, read_drift)) != 0) {
			window_err(err, yawfn); exit(err);
		}
#ifdef APP_DEBUG
	print_data_window(stderr, &coord);
#endif
	startt = detmax_startt (&coord, &rd, &pd, &yd);
	if (startt < 0.0) {
		fprintf (stderr, "%s: GPS and/or drift files not supplied!\n",
			 program_invocation_short_name);
		exit(EXIT_FAILURE);
	}
	imu[1]=0.0;
	while (imu[1] < startt) {
		if ((err = read_imufile(fpimu, imu)) != 1) {
			fprintf (stderr, "%s: first GPS time after IMU time\n",
				 program_invocation_short_name);
			exit(EXIT_FAILURE);
		}
	}
	err=0;
	while (err>=0 && read_imufile(fpimu, imu)!=EOF) {
		double rol = 0.0, pit = 0.0, yaw = 0.0;

		/* interpolate GPS to time */
		if (gpsfname != NULL) {
			err = update_data_window(&coord, imu[1], read_coord);
			if (err < 0) {
				window_err(err, coord.fname);
				exit(EXIT_FAILURE);
			}
#ifdef APP_DEBUG
			fprintf(stderr,"newtime %.6f\n", imu[1]);
			print_data_window(stderr, &coord);
#endif
			if (gps_interpol (&coord, imu[1], gpspos)!=0)
				exit(EXIT_FAILURE);
		}
		/* interpolate angle drift to time */
		if (rollfn != NULL) {
			err = update_data_window(&rd, imu[1], read_drift);
			if (err < 0) {
				window_err(err, rd.fname);
				exit(EXIT_FAILURE);
			}
			if (angle_interpol (&rd, imu[1], &rol)!=0) 
				exit(EXIT_FAILURE);
		}
		if (pitchfn != NULL) {
			err = update_data_window(&pd, imu[1], read_drift);
			if (err < 0) {
				window_err(err, pd.fname);
				exit(EXIT_FAILURE);
			}
			if (angle_interpol (&pd, imu[1], &pit)!=0) 
				exit(EXIT_FAILURE);
		}
		if (yawfn != NULL) {
			err = update_data_window(&yd, imu[1], read_drift);
			if (err < 0) {
				window_err(err, yd.fname);
				exit(EXIT_FAILURE);
			}
			if (angle_interpol (&yd, imu[1], &yaw)!=0) 
				exit(EXIT_FAILURE);
		}
		/* apply angle drift correction */
		angle_drift_corr (imu, rol, pit, yaw);
		/* apply offset coordinates w/ rotation */
		if (gpsfname != NULL) {
		  if (dooff) {
		    det_laser_pos (gpspos, antairframe, lasairframe,
				   imu[2]/RHO, imu[3]/RHO, imu[4]/RHO,
				   wgs84, laspos);
		  } else {
		    v3_copy (gpspos, laspos);
		  }
		}
		/* write output */
		fprintf (fpout, "%11.3f", imu[1]);
		if (gpsfname != NULL)
			fprintf (fpout, " %13.4f %13.4f %13.4f",
				 laspos[0], laspos[1], laspos[2]);
		fprintf (fpout, " %7.3f %7.3f %7.3f\n", imu[2], imu[3], imu[4]);
		
	}
	if (gpsfname != NULL) eliminate_data_window(&coord);
	if (rollfn != NULL) eliminate_data_window(&rd);
	if (pitchfn != NULL) eliminate_data_window(&pd);
	if (yawfn != NULL) eliminate_data_window(&yd);
}

/***************************************************************************/
int main (int argc, char **argv)
{
	double halfwin=1.0;

	extern char *program_invocation_name;
	extern char *program_invocation_short_name;
	char description[] = "merges GPS and iMAR IMU observations to nav file";
	char nonoptarg[] = " <imu_file> [<target>]\n   time types: gps 0, imu 4 drifts 4";
	FILE *fpimu, *fpout;
	int 
		ind_read,	/* index of first non option argument */
		quiet = FALSE,
		version_flag = FALSE;
	char *gpsfn = NULL;
	int applyoffset = TRUE;
	char *rollmed = NULL;
	char *pitchmed = NULL;
	char *yawmed = NULL;
	struct mainoptions myopts[] =
	{
		{0, 'V', NULL, 0, &version_flag, TRUE, "",
		 "version string."}
	      , {0, 'q', NULL, 0, &quiet, TRUE, "",
		 "run in quiet mode, e.g. do not comment output file."}
	      , {0, 'R', &rollmed, CMDL_CHAR, NULL, 0, " <roll_drift_file>",
		 "supply drift correction for roll, will be added to imu value."}
	      , {0, 'P', &pitchmed, CMDL_CHAR, NULL, 0, " <pitch_drift_file>",
		 "supply drift correction for pitch, will be added to imu value."}
	      , {0, 'Y', &yawmed, CMDL_CHAR, NULL, 0, " <yaw_drift_file>",
		 "supply drift correction for yaw, will be added to imu value."}
	      , {0, 'g', &gpsfn, CMDL_CHAR, NULL, 0, " <gps_file>",
		 "gps file to merge with INS data; if omitted, only drift will be added."}
	      , {0, 'L', NULL, 0, &applyoffset, FALSE, "",
		 "Do NOT apply laser offset; default is to apply laser offset."}
	      , {0, 'w', &halfwin, CMDL_DOUBLE, NULL, 0, " <half window size>",
 		 "half window size in seconds, default = 1.0sec"}
	};
	int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
      
	errinit (argv[0]);
	ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
					description, nonoptarg);
	if (version_flag) { 
		display_version();
		return EXIT_SUCCESS;
	}
	if (ind_read == -3)  return EXIT_SUCCESS;
	if (ind_read <0 || argc - ind_read < 1 ){
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
		return EXIT_FAILURE;
	}

/* where the work is done */
	if (!(fpimu=open_text_file(argv[ind_read],"r"))) return EXIT_FAILURE;
	if (argc - ind_read == 2) {
		fpout = open_text_file (argv[ind_read+1], "w");
		if (!fpout) return EXIT_FAILURE;
	} else {
		fpout = stdout;
	}

	if (quiet == FALSE) print_version (fpout, argv, argc);
	mknav1 (gpsfn, rollmed, pitchmed, yawmed,
		fpimu, fpout, applyoffset, halfwin, quiet);

	fclose (fpimu);
	fclose (fpout);
	return EXIT_SUCCESS;
}

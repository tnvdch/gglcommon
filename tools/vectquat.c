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

static void detangle (
	double d[3][3],
	double *roll,
	double *pitch,
	double *yaw)
{
	*pitch = asin (-d[0][2]);
	*yaw   = atan2 (d[0][1], d[0][0]);
	*roll  = atan2 (d[1][2], d[2][2]);
}


int vectquat (FILE *fpin, FILE *fpout, FILE *fperr, const int verbose)
{
	struct matrix *vg, *va, *vgvat, *ivgvgt, *rotmat, *rotmatt;
	double r, p, y, zeit, aq[4], qrm[3][3];
	int err, i, j;

	va    = m_init (3,4);
	vg    = m_init (3,4);
	vgvat = m_init (3,3);
	ivgvgt= m_init (3,3);
	rotmat= m_init (3,3);
	rotmatt=m_init (3,3);

  	vg->mat[0][0] =  3.921529; /* FRNT */
	vg->mat[1][0] =  0.296711; 
	vg->mat[2][0] = -0.127902;
  	vg->mat[0][1] = -5.351753; /* TAIL */
	vg->mat[1][1] = -0.233085;
	vg->mat[2][1] = -2.820545;
  	vg->mat[0][2] = -1.94856;  /* LEFT */
	vg->mat[1][2] =  3.4161;
	vg->mat[2][2] =  0.122967;
	vg->mat[0][3] =  2.84261;  /* RGHT */
	vg->mat[1][3] = -4.53794;
	vg->mat[2][3] =  2.82431;

	m_mult (vg, FALSE, vg, TRUE, ivgvgt);
        err = m_inv (ivgvgt);
	if (err != 0) return err;

	while((err=fscanf(fpin," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &zeit,
			  &(va->mat[0][0]), &(va->mat[1][0]), &(va->mat[2][0]),
			  &(va->mat[0][1]), &(va->mat[1][1]), &(va->mat[2][1]),
			  &(va->mat[0][2]), &(va->mat[1][2]), &(va->mat[2][2]),
			  &(va->mat[0][3]), &(va->mat[1][3]), &(va->mat[2][3]))
		) != EOF) {
		m_mult (vg, FALSE, va, TRUE, vgvat);
		m_mult (ivgvgt, FALSE, vgvat, FALSE, rotmatt);
/* 		for (i=0;i<3;i++) v3_norm (rotmat->mat[i],rotmat->mat[i]); */
		m_trans (rotmatt, rotmat);
		m_write (fperr, rotmat, 5,3,0,FALSE,TRUE);
		for (i=0;i<3;i++) for (j=0;j<3;j++) 
			qrm[i][j] = rotmat->mat[i][j];
		rotmat_to_quat (qrm, aq);
		quat_to_rotmat (aq, qrm);
		m3_write (fperr, qrm);
		detangle (qrm, &r, &p, &y);
                if (y<0.0) y+= 2.0*M_PI;
                fprintf (fpout, "Wc %12.8f %12.8f %12.8f %12.8f\n", zeit,
                         r/M_PI*180,p/M_PI*180,y/M_PI*180);
	}
	return 0;
}

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
	char description[] = "compute rotation matrix and angles from vectors";
	char nonoptarg[] = " <input_filename> [<input_filename ...]";
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
        if (isatty(STDIN_FILENO)) {
                fprintf (stderr, "You don't wanna type in the data, do you?\n\n");
                print_usage_info(FALSE, myopts, optcount, description, nonoptarg);
                return EXIT_FAILURE;
        }

	if (vectquat (stdin, stdout, stderr, verbose)!=0)
		fprintf (stderr, "error\n");

	return EXIT_SUCCESS;
}

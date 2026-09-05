#include <errno.h>      /*  the allround headers */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>      /*  character handling */
#include <string.h>     /*  string handling */
#include <fftw.h>
#include <math.h>
#include <common.h>     /*  our own institute's headers */

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"


/*********************************************************************/
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
                 "  maintainers: E. Favey <favey@geod.baug.ethz.ch>\n\n");
}

/*********************************************************************/

static void ff_error (int msgno)
{
	extern char *program_invocation_short_name;

	fprintf (stderr, "%s :", program_invocation_short_name);
	switch (msgno) {
	case 1:
		fprintf (stderr,"Error creating plan");
		break;
	case 2:
		fprintf (stderr,"Premature EOF or file read error: real part");
		break;
	case 3:
		fprintf (stderr,
			 "Premature EOF or file read error: imaginary part");
		break;
	default:
		fprintf (stderr,"ff error # %d", msgno);
	}
	fprintf (stderr, "\n");
	exit (EXIT_FAILURE);
}

/*********************************************************************/

static void call_fftw (
	FILE *fpin,
	FILE *fpout,
	const char *fname,
	const int numbering,
	const int n,
	const int backwards,
	const int complex,
	const int estimate)
{
	double reval, imval;
	int i;
	FFTW_COMPLEX *in = (FFTW_COMPLEX *) malloc (n* sizeof(FFTW_COMPLEX));
	FFTW_COMPLEX *out = (FFTW_COMPLEX *) malloc (n* sizeof(FFTW_COMPLEX));
	fftw_plan plan;
	double sqrtn = sqrt(n);

	plan = fftw_create_plan(n,
				((backwards) ? FFTW_BACKWARD : FFTW_FORWARD ),
				((estimate) ? FFTW_ESTIMATE : FFTW_MEASURE) );
	if (plan == NULL) ff_error (1);

	/* initialize array values */
	for (i = 0; i < n; ++i) {
		int err;
		reval = imval = 0.0;
		err = fscanf (fpin, " %lf", &reval);
		if (err != 1) ff_error (2);
		c_re(in[i]) = reval;
		if (complex) {
			err = fscanf (fpin, " %lf", &imval);
			if (err != 1) ff_error (3);
		}
		c_im(in[i]) = imval;
	}

	fftw (plan, 1, in, 1, 0, out, 1, 0);

	/* normalize */
	for (i = 0; i < n; ++i) {
		c_re(out[i]) = c_re(out[i])/sqrtn;
		c_im(out[i]) = c_im(out[i])/sqrtn;
	}

	/* output array */
	fprintf (fpout, "#  %s transform of %s\n", 
		 ((backwards)?"Backward":"Forward"), fname);
	fprintf (fpout, "# real part    imaginary   re^2+im^2\n");
	for (i = 0; i < n; ++i) {
	  if (numbering) fprintf (fpout, "%d ", i);
	  fprintf (fpout,"%12.5f %12.5f %12.5f\n", 
		   c_re(out[i]), c_im(out[i]), 
		   c_re(out[i])*c_re(out[i])+c_im(out[i])*c_im(out[i]));
	}
	fftw_destroy_plan (plan);
	free (in);
	free (out);
}

/*********************************************************************/

int main ( int argc, char **argv)
{
        int ind_read;           /* index of first non option argument */

	FILE *fp;
        int pr_version = FALSE;
	int n=0;
	int backwards = FALSE;
	int complex = FALSE;
	int estimate = FALSE;
	int linenr = FALSE;
	char *fname;
        struct mainoptions myopts[] =                   /* from cmdline.h */
        {
                /* a flag type entry */
                {0,'V', NULL, 0, &pr_version, TRUE, "",
                 "Version string."},
                {0,'n', &n, CMDL_INT, NULL, 0, " <number of values>",
                 "number of values for array."},
                {0,'b', NULL, 0, &backwards, TRUE, "",
                 "do backward fourier transform instead of forward."},
                {0,'c', NULL, 0, &complex, TRUE, "",
                 "read real and imaginary parts (default: only real)."},
                {0,'l', NULL, 0, &linenr, TRUE, "",
                 "number lines in first column (default: no numbers)."},
                {0,'E', NULL, 0, &estimate, TRUE, "",
                 "use estimate (default: measure) to create plan."}
        };
        int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
        char description[] = "Fast Fourier Transform 1 D";
        char nonoptarg[] = " [<input filename>]";
        extern char *program_invocation_name;           /* from errproc.h */
        extern char *program_invocation_short_name;     /* from errproc.h */

        errinit (argv[0]);                              /* from errproc.h */
        ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
                                        description, nonoptarg);
        if (pr_version) { display_version(); return EXIT_SUCCESS; }
        if (ind_read == -3) return EXIT_SUCCESS;
        if (ind_read < 0 || argc - ind_read > 1 || n < 1) {
                print_usage_info(FALSE, myopts, optcount, 
                                 description, nonoptarg);
                return EXIT_FAILURE;
        }
	if (argc - ind_read == 1) {
		fp = open_text_file (argv[ind_read], "r");
		if (!fp) return EXIT_FAILURE;
		fname = argv[ind_read];
	} else {
		fp = stdin;
		fname = strdup ("<STDIN>");
	}

	call_fftw (fp, stdout, fname, linenr, n, backwards, complex, estimate);

	if (fp) fclose (fp);
	return EXIT_SUCCESS;
}

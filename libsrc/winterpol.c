/*
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: winterpol.c,v 1.4 2000/02/25 12:57:58 ypeter Exp favey $
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 */
#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <time.h>	/*  for printing current time */
#include <math.h>	/*  mathematics */

#include "matcalc.h"
#include "window.h"
#include "winterpol.h"


/*******************************************************************/
/*
 *  RETURN VALUES: 0 on success, -1 on error
 */
int init_join (
	struct wfileinfo *fi,
	const char *fname)
{
	int ret;
	extern char *program_invocation_short_name;

	ret = init_data_window (fi->half_window_size,
				fi->blocksize,
				fi->num_cols,
				fname,
				&(fi->dw),
				fi->read_func);
	if (ret == EOF) {
		fprintf (stderr, "%s: %s: empty file\n",
			 program_invocation_short_name, fname);
		return -1;
	}
	if (ret != 0) {
		window_err(ret, fname);
		fprintf (stderr, "%s: %s: error %d\n",
			 program_invocation_short_name, fname, ret);
		return -1;
	}
#ifdef LIB_DEBUG
	print_data_window (stderr, &(fi->dw));
#endif
	return 0;
}


/*******************************************************************/
/*
 *  RETURN VALUES:	negative:	file update error
 *			positive:	no decent interpolation scheme
 *			zero:		updated and interpolated
 */
int interpol_join (
	struct wfileinfo *fi,	/* the file we operate on */
	const double thetime,	/* time of interpolation */
	const int num_coeff,	/* how many coefficients */
	const int donorm,	/* norm time between -1 .. +1 */
	const int cols,		/* column 1..cols are interpolated */
	double value[])		/* value[1] .. value[cols] hold results */
{
	int ret, i, brows;
	double ct;
	int xcol;

	fi->change = 0;
/*	do { */
	ret = update_data_window (&(fi->dw),
				  thetime,
				  fi->read_func);
	if (ret > 0) { fi->change += ret; fi->init = 0; }
/*	} while (ret > 0); */
	if (ret < 0) {
		window_err (ret, fi->dw.fname);
		return ret;
	}

 	if (donorm) {
		xcol = fi->dw.m->col-1;
		ct = ((thetime-fi->dw.tail_time)+(thetime-fi->dw.head_time)) /
			(fi->dw.head_time - fi->dw.tail_time);
	} else {
		xcol = 0;
		ct = thetime;
	}

	/* test if we have enough points before and after thetime */
	brows = 0;
	while (thetime > fi->dw.m->mat[brows][0] && brows<fi->dw.valid_count) 
		brows++;
	if (num_coeff > 2 && (brows < 2 || fi->dw.valid_count - brows < 2)) return 1;

	if (fi->change || !(fi->init)) {
		value[0] = ct;
 		if (donorm) for (i=0;i<fi->dw.valid_count;i++) {
 			fi->dw.m->mat[i][xcol] = 
 				((fi->dw.m->mat[i][0]-fi->dw.tail_time)+
				 (fi->dw.m->mat[i][0]-fi->dw.head_time) ) /
 				(fi->dw.head_time - fi->dw.tail_time);
 		}
		for (i=1;i<=cols;i++) {
			ret = m_polynom_approx (fi->dw.m,
						xcol, i, -1, -1, -1, 0, 
						fi->dw.valid_count-1,
						num_coeff, value[0], NULL, 
						fi->coeff[i], fi->qxx[i],
						&(fi->sigma[i]));
		}
		fi->init = 1;
	}

	for (i=1;i<=cols;i++) {
		value[i] = m_polynom_interpol (fi->coeff[i], num_coeff, ct-value[0]);
	}

	return 0;
}


/*******************************************************************/
/*
 *  main function joining and interpolating files

	/ * initialize * /
	if (outfname == NULL) {
		fpout = stdout;
	} else {
		fpout = open_text_file(outfname, "w");
	}
	if (!fpout) return -1;

	for (i=0;i<10;i++) {
		data->infile[0].coeff[i] = m_init (10,1);
		data->infile[1].coeff[i] = m_init (10,1);
		data->infile[0].qxx[i] = m_init (10,10);
		data->infile[1].qxx[i] = m_init (10,10);
	}
	if (init_join (&(data->infile[0]), gpsfname) != 0) return -1;
	if (init_join (&(data->infile[1]), plafname) != 0) return -1;
	if (init_join (&(data->infile[2]), lasfname) != 0) return -1;
	current_valgps = (double *) malloc (data->infile[0].num_cols*sizeof(double));
	current_valpla = (double *) malloc (data->infile[1].num_cols*sizeof(double));

	ttt = time(NULL);
	fprintf (fpout, "# %s: normalization = %s, Number of coeffs = %d, Date: %s",
		 program_invocation_short_name, (donorm)?"TRUE":"FALSE", num_coeff,
		 asctime(localtime(&ttt)));

	/ * working loop * /

	while (!feof(data->infile[2].dw.fpdata)) {
		double current_time;
		int fi1, fi2, k, drw;
		double phase;

		double tmpvec[3],  laserpos[3];
		double antairframe[3] = {-5.4097, 0.0574, -1.4996};
		double lasairframe[3] = {0.2802, 0.0040, -0.1200};
/ *  GPS Antenna in aircraft frame > -5.4097  0.0574 -1.4996  * /
/ *  Laser head in aircraft frame  >  0.2802  0.0040 -0.1200  * /

                / * read from high freq * /
		if (err != 1) {
			if (err != EOF)
				fprintf (stderr, "%s: error reading las file\n",
					 program_invocation_short_name);
			break;
		}
		current_time = data->infile[2].dw.m->mat[0][0];
		if ((fi1 = interpol_join (&(data->infile[0]),current_time,num_coeff,
					  donorm,3,current_valgps)) <0) {
			if (err != EOF)
				fprintf (stderr, "%s: error %d reading gps file\n",
					 program_invocation_short_name, fi1);
			break;
		}
		if (fi1 == 0 && fi2 ==0 && (fpout)) {
			fprintf (fpout, "1997 247 %15.8f %11.4f %11.4f %11.4f",
				 current_time, current_valgps[1],
				 current_valgps[2], current_valgps[3]);
			/ * calculate (r*e_L + b_A) * /
			data->infile[2].dw.m->mat[0][4] += ambig * 149.85;
			v3_multx(laserpos, data->infile[2].dw.m->mat[0][4], laserpos);
			v3_sub(lasairframe, antairframe, tmpvec);
			v3_add(laserpos, tmpvec, laserpos); 
			fprintf (fpout, " %13.7f %13.7f %13.7f %6.0f %2.0f",
				 laserpos[0], laserpos[1], laserpos[2],
				 data->infile[2].dw.m->mat[0][5],
				 data->infile[2].dw.m->mat[0][6]);
			if (data->verbose) {
				fprintf (stderr,
					 " sigmas %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n",
					 data->infile[0].sigma[1],
					 data->infile[0].sigma[2],
					 data->infile[0].sigma[3],
					 data->infile[1].sigma[1],
					 data->infile[1].sigma[2],
					 data->infile[1].sigma[3]);
			}
			fprintf (fpout, "\n");
		}

#ifdef LIB_DEBUG
		fprintf (stderr, "%d %d\n", fi1, fi2);
		print_data_window (stderr, &(data->infile[0].dw));
		print_data_window (stderr, &(data->infile[1].dw));
		print_data_window (stderr, &(data->infile[2].dw));
#endif

		err = 0;
	}

	/ * clean up * /
	eliminate_data_window (&(data->infile[0].dw));
	for (i=0;i<10;i++) {
		m_free(data->infile[0].coeff[i]);
		m_free(data->infile[1].coeff[i]);
		m_free(data->infile[0].qxx[i]);
		m_free(data->infile[1].qxx[i]);
	}
	fclose (fpout);

	return err;
}
*/

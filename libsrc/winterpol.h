/*
 *
 *  COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
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
 *  AUTHORS:  Etienne Favey <favey@geod.ethz.ch>
 *
 */
#ifndef GGL_WINTERPOL
#define GGL_WINTERPOL

#ifndef GGL_WINDOW
*** ERROR: include window.h prior to winterpol.h !
#endif

struct wfileinfo {
	FILE *fp;
	int blocksize;
	int num_cols;
	double half_window_size;
	struct data_window dw;
	int change;
	int init;
	struct matrix *coeff[10];
	struct matrix *qxx[10];
	double sigma[10];
	int (*read_func)(FILE *fpin, struct matrix *m, int row);
} ;



/*******************************************************************/
/*
 *  RETURN VALUES: 0 on success, -1 on error
 */
int init_join (
	struct wfileinfo *fi,
	const char *fname);



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
	double value[]);	/* value[1] .. value[cols] hold results */

#endif /* GGL_WINTERPOL */

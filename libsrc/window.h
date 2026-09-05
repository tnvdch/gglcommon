/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/window.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: file input and output functions including
 *            error handling
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: window.h,v 1.11 1999/09/08 09:12:18 favey Exp $
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
 *   HISTORY: 
 * $Log: window.h,v $
 * Revision 1.11  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.10  1999/06/30 12:17:30  favey
 * new function iszipped()
 *
 * Revision 1.9  1999/06/21 12:22:44  favey
 * bugfix in get_fname_wopath
 *
 * Revision 1.8  1999/06/21 09:47:52  favey
 * added file name to struct for easier logging
 *
 * Revision 1.7  1999/06/03 10:09:13  cocard
 * print_data_window: 2 format parameters added
 *
 * Revision 1.6  1999/04/08 13:32:52  favey
 * bugfix in time_index_data_window
 *
 * Revision 1.5  1998/10/06 16:19:02  cocard
 *  function time_index_data_window added
 *
 * Revision 1.4  1998/10/02 18:47:40  cocard
 * major changes to readfunction
 *
 * Revision 1.3  1998/07/21 12:40:31  favey
 * added error function
 *
 * Revision 1.2  1998/05/04 07:02:34  favey
 * *** empty log message ***
 *
 * Revision 1.1  1998/04/29 10:29:13  favey
 * Initial revision
 *
 *
 */

#ifndef GGL_MATCALC
*** ERROR: include matcalc.h prior to window.h !
#endif

#ifndef GGL_WINDOW
#define GGL_WINDOW


/*------------------------------------------------------------------
 *  Structure to hold file pointer, matrix pointer, time values, etc.
 */ 

struct data_window {
	double mid_time;	/* mid time of actual window */
	double tail_time;	/* start time of actual window */
	double head_time;	/* end time of actual window */
	double halfwinsize;	/* deltatime = (endtime - starttime)/2 */
	int bufsize;		/* maximal number of rows for a window */
	int valid_count;	/* actual number of rows */
	FILE *fpdata;		/* file pointer to data file */
	char *fname;		/* file name of data file */
	int ispipe;		/* flag whether fpdata is a pipe */
	struct matrix *m;	/* struct containing the data */
};

/*------------------------------------------------------------------
 *  Declaration of a reading function to be called
 *
 *  This function has to be written by the user. It must fill the 
 * specified row of the matrix. The value in the first column is 
 * considered to be the time.
 *
 *  The return value should be one of the following:
 *       1         when the reading of one line was successful
 *       0         when a file input/output error occurred
 *     EOF         when we reached the end of file
 *
 *  here an example:
 *    int read_data(FILE *fpin, struct matrix *m, int row)
 */


/*------------------------------------------------------------------
 *  initialize windowing structure
 *
 *  The return values may be:
 *       = 0       when the update was successful
 *       < 0       error: see function window_err
 */
int init_data_window (
	const double half_window_size,
	const int buffer_size,
	const int numcols,
	const char *datafname,
	struct data_window *w,
	int (*read_new_element) (FILE *fpin, struct matrix *m, int row));

/*------------------------------------------------------------------
 *  update windowing structure
 * 
 *  Moves window forward over file by at most buffersize values
 *  or until tail and head times 
 *
 *  The return values may be:
 *      >= 0       the amount of lines read
 *       < 0       error: see function window_err
 */
int update_data_window (
	struct data_window *w,
	const double new_time,
	int (*read_new_element)(FILE *fpin, struct matrix *m, int row));

/*------------------------------------------------------------------
 *  eliminate windowing structure
 *  Frees the matrix and close the file, sets all values to zero.
 */
void eliminate_data_window (struct data_window *w);

/*------------------------------------------------------------------
 *  print structure window (for debugging purposes)
 */
void print_data_window (FILE *fpout, struct data_window *w, int fmt1, int fmt2);

/*------------------------------------------------------------------
 * returns index i of row for  t(i) <= t < t(i+1)
 */
int time_index_data_window(struct data_window *w, double t);

/*------------------------------------------------------------------
 *  print appropriate error message
 *
 *    -1: end of file
 *    -2: I/O error or incorrect matrix update by your function
 *    -3: necessary memory could not be allocated
 *    -4: time not in window
 *    -5: buffer too small
 */
void window_err (const int errval, const char *fname);

#endif
